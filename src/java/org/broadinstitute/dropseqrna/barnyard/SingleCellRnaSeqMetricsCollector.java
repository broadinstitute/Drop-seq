/*
 * MIT License
 *
 * Copyright 2017 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package org.broadinstitute.dropseqrna.barnyard;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.broadinstitute.dropseqrna.TranscriptomeException;
import org.broadinstitute.dropseqrna.annotation.GeneAnnotationReader;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;
import org.broadinstitute.dropseqrna.utils.StringTagComparator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamRecordSortingIteratorFactory;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.samtools.util.ProgressLogger;
import picard.analysis.MetricAccumulationLevel;
import picard.analysis.RnaSeqMetrics;
import picard.analysis.directed.RnaSeqMetricsCollector;
import picard.annotation.Gene;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.metrics.PerUnitMetricCollector;

/**
 * An adaptation of the Picard RnaSeqMetricsCollector to collect per-cell data.  In particular, the exon/intron/genic/intragenic/rRNA levels.
 * @author nemesh
 *
 */

@CommandLineProgramProperties(
        usage = "An adaptation of the Picard RnaSeqMetricsCollector to collect per-cell data.  In particular, the exon/intron/genic/intragenic/rRNA levels" +
        		" This program looks at the mapping from each of the reads in both genomic and library space, and selects the better mapping.",
        usageShort = "Measures the intron/exon/genic/intergenic/rRNA levels of each cell.",
        programGroup = DropSeq.class
)
public class SingleCellRnaSeqMetricsCollector extends CommandLineProgram {

	private static final Log log = Log.getInstance(SingleCellRnaSeqMetricsCollector.class);

	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.")
	public File INPUT;

	@Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file of per-cell exonic/intronic/genic/intergenic/rRNA levels.  This supports zipped formats like gz and bz2.")
	public File OUTPUT;

	@Option(doc="The cell barcode tag.  If there are no reads with this tag, the program will assume that all reads belong to the same cell and process in single sample mode.")
	public String CELL_BARCODE_TAG="XC";

	@Option(doc="Gene annotations in refFlat or GTF format.")
	public File ANNOTATIONS_FILE;

    @Option(doc="Location of rRNA sequences in genome, in interval_list format.  " +
            "If not specified no bases will be identified as being ribosomal.  " +
            "Format described here: http://picard.sourceforge.net/javadoc/net/sf/picard/util/IntervalList.html", optional = true)
    public File RIBOSOMAL_INTERVALS;

    // for backwards compatability, if the strand isn't set, set it to none.
    // TODO should this default be set to FIRST_READ?
    @Option(shortName = "STRAND", doc="For strand-specific library prep. " +
            "For unpaired reads, use FIRST_READ_TRANSCRIPTION_STRAND if the reads are expected to be on the transcription strand.")
    public RnaSeqMetricsCollector.StrandSpecificity STRAND_SPECIFICITY = RnaSeqMetricsCollector.StrandSpecificity.NONE;

    @Option(doc="This percentage of the length of a fragment must overlap one of the ribosomal intervals for a read or read pair by this must in order to be considered rRNA.")
    public double RRNA_FRAGMENT_PERCENTAGE = 0.8;

    @Option(doc="Number of cells that you think are in the library. The top NUM_CORE_BARCODES will be reported in the output.", mutex={"CELL_BC_FILE"})
	public Integer NUM_CORE_BARCODES=null;

    @Option(doc="Override NUM_CORE_BARCODES, and process reads that have the cell barcodes in this file instead.  When supplied, output is ordered to match the input barcode ordering. The file has 1 column with no header.", mutex={"NUM_CORE_BARCODES"})
	public File CELL_BC_FILE=null;

    @Option(doc="The map quality of the read to be included for determining which cells will be measured.")
	public Integer READ_MQ=10;

    @Option(doc="If specified, count bases that align to this sequence separately from other categories")
    public List<String> MT_SEQUENCE;

    @Override
	protected int doWork() {
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);

		for (final String mtSequence : MT_SEQUENCE) {
			final SAMSequenceRecord samSequenceRecord =
					SamReaderFactory.makeDefault().open(INPUT).getFileHeader().getSequence(mtSequence);
			if (samSequenceRecord == null)
				throw new RuntimeException("MT_SEQUENCE '" + mtSequence + "' is not found in sequence dictionary in " + INPUT.getAbsolutePath());
		}

		List<String> cellBarcodes = getCellBarcodes(this.CELL_BC_FILE, this.INPUT, this.CELL_BARCODE_TAG, this.READ_MQ, this.NUM_CORE_BARCODES);
		RnaSeqMetricsCollector collector = getRNASeqMetricsCollector(this.CELL_BARCODE_TAG, cellBarcodes, this.INPUT, this.STRAND_SPECIFICITY, this.RRNA_FRAGMENT_PERCENTAGE, this.READ_MQ, this.ANNOTATIONS_FILE, this.RIBOSOMAL_INTERVALS);
		final MetricsFile<RnaSeqMetrics, Integer> file = getMetricsFile();
		log.info("Adding metrics to file.  This may take a while, with no progress messages.");
    	collector.addAllLevelsToFile(file);

    	BufferedWriter b = IOUtil.openFileForBufferedWriting(OUTPUT);
    	file.write(b);
    	try {
			b.close();
		} catch (IOException io) {
			throw new TranscriptomeException("Problem writing file", io);
		}
    	return 0;
    }

    /**
     * If there's a cell barcode file that is non-null, use that to get a list of cell barcodes.
     * Otherwise, gather up the top <numCoreBarcodes> cells ordered by number of reads.
     */
    private List<String> getCellBarcodes(final File cellBCFile, final File bamFile, final String cellBarcodeTag, final int readMQ, final Integer numCoreBarcodes) {
    	if (cellBCFile!=null) {
    		List<String>cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(cellBCFile);
    		log.info("Found " + cellBarcodes.size()+ " cell barcodes in file");
    		return cellBarcodes;
    	}
    	BarcodeListRetrieval u = new BarcodeListRetrieval();
        return u.getListCellBarcodesByReadCount (bamFile, cellBarcodeTag, readMQ, null, numCoreBarcodes);

    }

    RnaSeqMetricsCollector getRNASeqMetricsCollector(final String cellBarcodeTag, final List<String> cellBarcodes, final File inBAM,
    		final RnaSeqMetricsCollector.StrandSpecificity strand, final double rRNAFragmentPCT, final int readMQ,
    		final File annotationsFile, final File rRNAIntervalsFile) {

    	CollectorFactory factory = new CollectorFactory(inBAM, strand, rRNAFragmentPCT, annotationsFile, rRNAIntervalsFile);
		RnaSeqMetricsCollector collector=  factory.getCollector(cellBarcodes);
		List<SAMReadGroupRecord> rg = factory.getReadGroups(cellBarcodes);

        // iterate by cell barcodes.  Skip all the reads without cell barcodes.
		CloseableIterator<SAMRecord> iter = getReadsInTagOrder (inBAM, cellBarcodeTag, rg, cellBarcodes, readMQ);

        ProgressLogger p = new ProgressLogger(log, 1000000, "Accumulating metrics");
		while (iter.hasNext()) {
			SAMRecord r = iter.next();
			String cellBarcode = r.getStringAttribute(cellBarcodeTag);
			r.setAttribute("RG", cellBarcode);
            p.record(r);
	    	collector.acceptRecord(r, null);
		}

		collector.finish();
		return (collector);
    }


    /**
     * Sets up the reads in cell barcode order.
     * Only adds reads that pass the map quality and are in the set of cell barcodes requested.
     *
     * I've tried adapting this to the TagOrderIterator API, but it seems like I need to add the read groups to the header of the temporary BAM that gets
     * iterated on or this doesn't work.
     */
    private CloseableIterator<SAMRecord> getReadsInTagOrder (final File bamFile, final String primaryTag,
                                                             final List<SAMReadGroupRecord> rg,
                                                             final List<String> allCellBarcodes, final int mapQuality) {

		SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
		SAMSequenceDictionary dict= reader.getFileHeader().getSequenceDictionary();
		List<SAMProgramRecord> programs =reader.getFileHeader().getProgramRecords();

		final Set<String> cellBarcodeSet = new HashSet<> (allCellBarcodes);

		final SAMFileHeader writerHeader = new SAMFileHeader();
		// reader.getFileHeader().setReadGroups(rg);
		for (SAMReadGroupRecord z: rg) {
			reader.getFileHeader().addReadGroup(z);
			writerHeader.addReadGroup(z);
		}
        writerHeader.setSortOrder(SAMFileHeader.SortOrder.queryname);
        writerHeader.setSequenceDictionary(dict);
        for (SAMProgramRecord spr : programs)
			writerHeader.addProgramRecord(spr);

        // This not only filters, but sets the RG attribute on reads it allows through.
        final FilteredIterator<SAMRecord> rgAddingFilter = new FilteredIterator<SAMRecord>(reader.iterator()) {
            @Override
            public boolean filterOut(final SAMRecord r) {
                String cellBarcode = r.getStringAttribute(primaryTag);
                if (cellBarcodeSet.contains(cellBarcode) & r.getMappingQuality() >= mapQuality) {
                    r.setAttribute("RG", cellBarcode);
                    return false;
                } else
					return true;
            }
        };

        ProgressLogger p = new ProgressLogger(log, 1000000, "Preparing reads in core barcodes");
        CloseableIterator<SAMRecord> sortedIterator = SamRecordSortingIteratorFactory.create(writerHeader, rgAddingFilter, new StringTagComparator(primaryTag), p);

		log.info("Sorting finished.");
		return (sortedIterator);
	}

    private class CollectorFactory {
    	final OverlapDetector<Gene> geneOverlapDetector;
    	final Long ribosomalBasesInitialValue;
    	final OverlapDetector<Interval> ribosomalSequenceOverlapDetector;
    	final HashSet<Integer> ignoredSequenceIndices;
    	final RnaSeqMetricsCollector.StrandSpecificity specificity;
    	final double rnaFragPct;

    	public CollectorFactory (final File bamFile, final RnaSeqMetricsCollector.StrandSpecificity specificity, final double rnaFragPct, final File annotationsFile, final File ribosomalIntervals) {
    		this.specificity=specificity;
    		this.rnaFragPct=rnaFragPct;
    		SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
    		geneOverlapDetector = GeneAnnotationReader.loadAnnotationsFile(annotationsFile, reader.getFileHeader().getSequenceDictionary());
            log.info("Loaded " + geneOverlapDetector.getAll().size() + " genes.");
            ribosomalBasesInitialValue = ribosomalIntervals != null ? 0L : null;
            ribosomalSequenceOverlapDetector = RnaSeqMetricsCollector.makeOverlapDetector(bamFile, reader.getFileHeader(), ribosomalIntervals, log);
            ignoredSequenceIndices = RnaSeqMetricsCollector.makeIgnoredSequenceIndicesSet(reader.getFileHeader(), new HashSet<String>());
            CloserUtil.close(reader);
    	}

    	public RnaSeqMetricsCollector getCollector(final List<String> cellBarcodes) {
    		List<SAMReadGroupRecord> readGroups =  getReadGroups(cellBarcodes);
    		return new RnaSeqMtMetricsCollector(CollectionUtil.makeSet(MetricAccumulationLevel.READ_GROUP), readGroups,
                    ribosomalBasesInitialValue, geneOverlapDetector, ribosomalSequenceOverlapDetector,
                    ignoredSequenceIndices, 500, specificity, this.rnaFragPct, false);
    	}

    	public List<SAMReadGroupRecord> getReadGroups(final List<String> cellBarcodes) {
    		List<SAMReadGroupRecord> g = new ArrayList<>(cellBarcodes.size());
    		for (String id: cellBarcodes) {
    			SAMReadGroupRecord rg = new SAMReadGroupRecord(id);
    			rg.setLibrary(id);
    		    rg.setPlatform(id);
    		    rg.setSample(id);
    		    rg.setPlatformUnit(id);
    			g.add(rg);
    		}
    		return (g);


    	}
    }

    public static class RnaSeqMtMetrics
            extends RnaSeqMetrics {
        public long MT_BASES;
        public double PCT_MT_BASES;
    }

	private class RnaSeqMtMetricsCollector extends RnaSeqMetricsCollector {
        public RnaSeqMtMetricsCollector(final Set<MetricAccumulationLevel> accumulationLevels,
                                        final List<SAMReadGroupRecord> samRgRecords,
                                        final Long ribosomalBasesInitialValue,
                                        final OverlapDetector<Gene> geneOverlapDetector,
                                        final OverlapDetector<Interval> ribosomalSequenceOverlapDetector,
                                        final HashSet<Integer> ignoredSequenceIndices, final int minimumLength,
                                        final StrandSpecificity strandSpecificity,
                                        final double rrnaFragmentPercentage,
                                        final boolean collectCoverageStatistics) {
            super(accumulationLevels, samRgRecords, ribosomalBasesInitialValue, geneOverlapDetector,
                    ribosomalSequenceOverlapDetector, ignoredSequenceIndices, minimumLength, strandSpecificity,
                    rrnaFragmentPercentage, collectCoverageStatistics);
        }

        @Override
        protected PerUnitMetricCollector<RnaSeqMetrics, Integer, SAMRecord> makeChildCollector(final String sample, final String library, final String readGroup) {
            // TODO: Fix last argument
            return new PerUnitRnaSeqMtMetricsCollector(sample, library, readGroup, this.ribosomalInitialValue);
        }

        private class PerUnitRnaSeqMtMetricsCollector extends PerUnitRnaSeqMetricsCollector {

            public PerUnitRnaSeqMtMetricsCollector(final String sample, final String library, final String readGroup, final Long ribosomalBasesInitialValue) {
                super(new RnaSeqMtMetrics(), sample, library, readGroup, ribosomalBasesInitialValue);
            }

			private RnaSeqMtMetrics castMetrics() {
				return (RnaSeqMtMetrics)metrics;
			}

            @Override
            public void acceptRecord(final SAMRecord rec) {
                if (MT_SEQUENCE.contains(rec.getReferenceName()) && !rec.getReadFailsVendorQualityCheckFlag() &&
                        !rec.isSecondaryOrSupplementary() && !rec.getReadUnmappedFlag()) {

                    metrics.PF_BASES += rec.getReadLength();
                    final int numAlignedBases = getNumAlignedBases(rec);
                    castMetrics().MT_BASES += numAlignedBases;
                    metrics.PF_ALIGNED_BASES += numAlignedBases;
                } else
					super.acceptRecord(rec);
            }

            @Override
            public void finish() {
                super.finish();
                if (metrics.PF_ALIGNED_BASES > 0)
					castMetrics().PCT_MT_BASES = castMetrics().MT_BASES / (double) metrics.PF_ALIGNED_BASES;
            }

        }
    }

    /** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new SingleCellRnaSeqMetricsCollector().instanceMain(args));
	}
}
