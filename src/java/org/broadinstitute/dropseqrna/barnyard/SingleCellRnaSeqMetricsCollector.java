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
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.dropseqrna.TranscriptomeException;
import org.broadinstitute.dropseqrna.annotation.GeneAnnotationReader;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
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
import picard.cmdline.StandardOptionDefinitions;
import picard.metrics.PerUnitMetricCollector;

/**
 * An adaptation of the Picard RnaSeqMetricsCollector to collect per-cell data.  In particular, the exon/intron/genic/intragenic/rRNA levels.
 * @author nemesh
 *
 */

@CommandLineProgramProperties(
        summary = "An adaptation of the Picard RnaSeqMetricsCollector to collect per-cell data.  In particular, the exon/intron/genic/intragenic/rRNA levels" +
        		" This program looks at the mapping from each of the reads in both genomic and library space, and selects the better mapping.",
        oneLineSummary = "Measures the intron/exon/genic/intergenic/rRNA levels of each cell.",
        programGroup = DropSeq.class
)
public class SingleCellRnaSeqMetricsCollector extends CommandLineProgram {

    private static final Log log = Log.getInstance(SingleCellRnaSeqMetricsCollector.class);
    private static final int MINIMUM_TRANSCRIPT_LENGTH = 500;

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze. This argument can accept wildcards, or a file with the suffix .bam_list that contains the locations of multiple BAM files", minElements = 1)
	public List<File> INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file of per-cell exonic/intronic/genic/intergenic/rRNA levels.  This supports zipped formats like gz and bz2.")
	public File OUTPUT;

	@Argument(doc="The cell barcode tag.  If there are no reads with this tag, the program will assume that all reads belong to the same cell and process in single sample mode.")
	public String CELL_BARCODE_TAG="XC";

	@Argument(doc="Gene annotations in refFlat or GTF format.")
	public File ANNOTATIONS_FILE;

    @Argument(doc="Location of rRNA sequences in genome, in interval_list format.  " +
            "If not specified no bases will be identified as being ribosomal.  " +
            "Format described here: http://picard.sourceforge.net/javadoc/net/sf/picard/util/IntervalList.html", optional = true)
    public File RIBOSOMAL_INTERVALS;

    // For backwards compatibility the default is NONE (no strand-specific metrics), not FIRST_READ_TRANSCRIPTION_STRAND.
    @Argument(shortName = "STRAND", doc="For strand-specific library prep. " +
            "For unpaired reads, use FIRST_READ_TRANSCRIPTION_STRAND if the reads are expected to be on the transcription strand.")
    public RnaSeqMetricsCollector.StrandSpecificity STRAND_SPECIFICITY = RnaSeqMetricsCollector.StrandSpecificity.NONE;

    @Argument(doc="This percentage of the length of a fragment must overlap one of the ribosomal intervals for a read or read pair by this must in order to be considered rRNA.")
    public double RRNA_FRAGMENT_PERCENTAGE = 0.8;

    @Argument(doc="Number of cells that you think are in the library. The top NUM_CORE_BARCODES will be reported in the output.", mutex={"CELL_BC_FILE"})
	public Integer NUM_CORE_BARCODES=null;

    @Argument(doc="Override NUM_CORE_BARCODES, and process reads that have the cell barcodes in this file instead.  When supplied, output is ordered to match the input barcode ordering. The file has 1 column with no header.", mutex={"NUM_CORE_BARCODES"})
	public File CELL_BC_FILE=null;

    @Argument(doc="The map quality of the read to be included for determining which cells will be measured.")
	public Integer READ_MQ=10;

    @Argument(doc="If specified, count bases that align to this sequence separately from other categories",
            optional = true)
    public List<String> MT_SEQUENCE;

    @Argument(doc="Number of threads.  With 1 thread, reads are decoded and accumulated on the calling thread.  With more, " +
            "the calling thread reads, decodes and filters the input BAM(s) and the other threads accumulate the per-cell metrics.")
    public int NUM_THREADS = 1;

    @Hidden
    @Argument(doc="Threads used to read and decode the input.  Must be given together with WORKER_THREADS, which then replace NUM_THREADS.", optional = true)
    public Integer READER_THREADS = null;

    @Hidden
    @Argument(doc="Threads used to accumulate metrics (0 = on the calling thread).  Must be given together with READER_THREADS.", optional = true)
    public Integer WORKER_THREADS = null;

	@Override
	protected int doWork() {
		validateAndExpandInputs();

		List<String> cellBarcodes = getCellBarcodes(this.CELL_BC_FILE, this.INPUT, this.CELL_BARCODE_TAG, this.READ_MQ, this.NUM_CORE_BARCODES);
		final int[] threads = getReaderAndWorkerThreads();
		log.info("Reader threads: " + threads[0] + ", worker threads: " + threads[1]);
		RnaSeqMetricsCollector collector = getRNASeqMetricsCollector(this.CELL_BARCODE_TAG, cellBarcodes, this.INPUT, this.STRAND_SPECIFICITY, this.RRNA_FRAGMENT_PERCENTAGE, this.READ_MQ, this.ANNOTATIONS_FILE, this.RIBOSOMAL_INTERVALS, threads[0], threads[1]);
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
    List<String> getCellBarcodes(final File cellBCFile, final List<File> bamFiles, final String cellBarcodeTag, final int readMQ, final Integer numCoreBarcodes) {
    	if (cellBCFile!=null) {
    		List<String>cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(cellBCFile);
    		log.info("Found " + cellBarcodes.size()+ " cell barcodes in file");
    		return cellBarcodes;
    	}
    	BarcodeListRetrieval u = new BarcodeListRetrieval();
        return u.getListCellBarcodesByReadCount (bamFiles, cellBarcodeTag, readMQ, null, numCoreBarcodes);

    }

    RnaSeqMetricsCollector getRNASeqMetricsCollector(final String cellBarcodeTag, final List<String> cellBarcodes, final List<File> inBAMs,
    		final RnaSeqMetricsCollector.StrandSpecificity strand, final double rRNAFragmentPCT, final int readMQ,
    		final File annotationsFile, final File rRNAIntervalsFile) {
    	return getRNASeqMetricsCollector(cellBarcodeTag, cellBarcodes, inBAMs, strand, rRNAFragmentPCT, readMQ, annotationsFile, rRNAIntervalsFile, 1, 0);
    }

    /**
     * Reads are not sorted: each read goes to the collector of its own cell, and the per-cell metrics do not depend on
     * the order reads arrive.  Only reads that pass the map quality and are in the set of cell barcodes requested are used.
     * @param readerThreads threads for reading/decoding the input
     * @param workerThreads threads accumulating metrics.  0 means on the calling thread.  Each cell is always handled by one thread.
     */
    RnaSeqMetricsCollector getRNASeqMetricsCollector(final String cellBarcodeTag, final List<String> cellBarcodes, final List<File> inBAMs,
    		final RnaSeqMetricsCollector.StrandSpecificity strand, final double rRNAFragmentPCT, final int readMQ,
    		final File annotationsFile, final File rRNAIntervalsFile, final int readerThreads, final int workerThreads) {

    	CollectorFactory factory = new CollectorFactory(inBAMs.get(0), strand, rRNAFragmentPCT, annotationsFile, rRNAIntervalsFile);
		final RnaSeqMtMetricsCollector collector = factory.getCollector(cellBarcodes);
		// Reads go straight to the collector of their cell, so no read group is needed on the reads or in the header.
		final SingleCellRnaSeqMetricsAccumulator accumulator = new SingleCellRnaSeqMetricsAccumulator(collector.getCollectorsByCell(), workerThreads);

		final SamHeaderAndIterator headerAndIter = SamFileMergeUtil.mergeInputPaths(IOUtil.filesToPaths(inBAMs), false, SamReaderFactory.makeDefault(), readerThreads);
		final CloseableIterator<SAMRecord> iter = headerAndIter.iterator;
		final ProgressLogger p = new ProgressLogger(log, 1000000, "Accumulating metrics");
		while (iter.hasNext()) {
			final SAMRecord r = iter.next();
			if (r.getMappingQuality() < readMQ)
				continue;
			final String cellBarcode = r.getStringAttribute(cellBarcodeTag);
			if (cellBarcode != null && accumulator.accept(cellBarcode, r))
				p.record(r);
		}
		accumulator.finish();
		CloserUtil.close(headerAndIter.iterator);

		collector.finish();
		return (collector);
    }

    private class CollectorFactory {
    	final OverlapDetector<Gene> geneOverlapDetector;
    	final Set<Gene> genesWithLongEnoughTranscripts;
    	final Long ribosomalBasesInitialValue;
    	final OverlapDetector<Interval> ribosomalSequenceOverlapDetector;
    	final HashSet<Integer> ignoredSequenceIndices;
    	final RnaSeqMetricsCollector.StrandSpecificity specificity;
    	final double rnaFragPct;

    	public CollectorFactory (final File bamFile, final RnaSeqMetricsCollector.StrandSpecificity specificity, final double rnaFragPct, final File annotationsFile, final File ribosomalIntervals) {
    		this.specificity=specificity;
    		this.rnaFragPct=rnaFragPct;
    		SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
    		geneOverlapDetector = GeneAnnotationReader.loadAnnotationsFile(annotationsFile,
					reader.getFileHeader().getSequenceDictionary(), VALIDATION_STRINGENCY);
    		// This is a time-consuming call, so invoke it once here so that it isn't invoke over and over by RnaSeqMetricsCollector
    		genesWithLongEnoughTranscripts = geneOverlapDetector.getAll();
            log.info("Loaded " + genesWithLongEnoughTranscripts.size() + " genes.");
            ribosomalBasesInitialValue = ribosomalIntervals != null ? 0L : null;
            ribosomalSequenceOverlapDetector = RnaSeqMetricsCollector.makeOverlapDetector(bamFile, reader.getFileHeader(), ribosomalIntervals, log);
            ignoredSequenceIndices = RnaSeqMetricsCollector.makeIgnoredSequenceIndicesSet(reader.getFileHeader(), new HashSet<String>());
            CloserUtil.close(reader);
    	}

    	public RnaSeqMtMetricsCollector getCollector(final List<String> cellBarcodes) {
    		List<SAMReadGroupRecord> readGroups =  getReadGroups(cellBarcodes);
    		return new RnaSeqMtMetricsCollector(CollectionUtil.makeSet(MetricAccumulationLevel.READ_GROUP), readGroups,
                    ribosomalBasesInitialValue, geneOverlapDetector, ribosomalSequenceOverlapDetector,
                    ignoredSequenceIndices, MINIMUM_TRANSCRIPT_LENGTH, specificity, this.rnaFragPct, false,
                    genesWithLongEnoughTranscripts);
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

	private class RnaSeqMtMetricsCollector extends RnaSeqMetricsCollector {
        private Set<Gene> genesWithLongEnoughTranscripts;
        // Deliberately has no initializer: makeChildCollector is called from the superclass constructor, before
        // field initializers run, and an initializer would discard the children already created.
        private Map<String, PerUnitMetricCollector<RnaSeqMetrics, Integer, SAMRecord>> collectorsByCell;

        public RnaSeqMtMetricsCollector(final Set<MetricAccumulationLevel> accumulationLevels,
                                        final List<SAMReadGroupRecord> samRgRecords,
                                        final Long ribosomalBasesInitialValue,
                                        final OverlapDetector<Gene> geneOverlapDetector,
                                        final OverlapDetector<Interval> ribosomalSequenceOverlapDetector,
                                        final HashSet<Integer> ignoredSequenceIndices, final int minimumLength,
                                        final StrandSpecificity strandSpecificity,
                                        final double rrnaFragmentPercentage,
                                        final boolean collectCoverageStatistics,
                                        final Set<Gene> genesWithLongEnoughTranscripts) {
            super(accumulationLevels, samRgRecords, ribosomalBasesInitialValue, geneOverlapDetector,
                    ribosomalSequenceOverlapDetector, ignoredSequenceIndices, minimumLength, strandSpecificity,
                    rrnaFragmentPercentage, collectCoverageStatistics);
            this.genesWithLongEnoughTranscripts = genesWithLongEnoughTranscripts;
        }

        @Override
        protected PerUnitMetricCollector<RnaSeqMetrics, Integer, SAMRecord> makeChildCollector(final String sample, final String library, final String readGroup) {
            // Same as Picard's own makeChildCollector: ribosomalInitialValue is set before the superclass calls this.
            final PerUnitMetricCollector<RnaSeqMetrics, Integer, SAMRecord> child =
                    new PerUnitRnaSeqMtMetricsCollector(sample, library, readGroup, this.ribosomalInitialValue);
            if (collectorsByCell == null)
                collectorsByCell = new HashMap<>();
            collectorsByCell.put(readGroup, child);
            return child;
        }

        /** The per-cell collectors (read group = cell barcode), so reads can be given to them directly. */
        Map<String, PerUnitMetricCollector<RnaSeqMetrics, Integer, SAMRecord>> getCollectorsByCell() {
            return collectorsByCell;
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
                if (MT_SEQUENCE!=null &&  MT_SEQUENCE.contains(rec.getReferenceName()) && !rec.getReadFailsVendorQualityCheckFlag() &&
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
                // PCT_RIBOSOMAL_BASES should never be null.
                if (metrics.PCT_RIBOSOMAL_BASES==null) metrics.PCT_RIBOSOMAL_BASES=0d;
                if (metrics.PF_ALIGNED_BASES > 0)
					castMetrics().PCT_MT_BASES = castMetrics().MT_BASES / (double) metrics.PF_ALIGNED_BASES;
            }

            @Override
            protected Set<Gene> getGenesForPickTranscripts() {
                return genesWithLongEnoughTranscripts;
            }
        }
    }

    @Override
    protected String[] customCommandLineValidation() {
        if (NUM_THREADS < 1)
            return new String[] { "NUM_THREADS must be >= 1" };
        if ((READER_THREADS == null) != (WORKER_THREADS == null))
            return new String[] { "READER_THREADS and WORKER_THREADS must be given together" };
        if (READER_THREADS != null && (READER_THREADS < 1 || WORKER_THREADS < 0))
            return new String[] { "READER_THREADS must be >= 1 and WORKER_THREADS >= 0" };
        return super.customCommandLineValidation();
    }

    /**
     * @return {reader threads, worker threads}.  Zero workers means metrics are accumulated on the calling thread.
     */
    int[] getReaderAndWorkerThreads() {
        if (READER_THREADS != null)
            return new int[] { READER_THREADS, WORKER_THREADS };
        // The accumulating workers are the bottleneck, and a second decode thread did not pay for itself.
        return new int[] { 1, NUM_THREADS - 1 };
    }

    /**
     * Expands wildcards / .bam_list inputs, checks that input and output files are usable, and checks that
     * any MT_SEQUENCE is in the sequence dictionary of the first input.
     */
    private void validateAndExpandInputs() {
		this.INPUT = FileListParsingUtils.expandFileList(INPUT);
		for (final File f : INPUT) IOUtil.assertFileIsReadable(f);
		IOUtil.assertFileIsWritable(OUTPUT);
		if (RIBOSOMAL_INTERVALS!=null) IOUtil.assertFileIsReadable(RIBOSOMAL_INTERVALS);

        if (MT_SEQUENCE == null) {
            MT_SEQUENCE = new ArrayList<>();
        }
		for (final String mtSequence : MT_SEQUENCE) {
			final SAMSequenceRecord samSequenceRecord =
					SamReaderFactory.makeDefault().open(INPUT.get(0)).getFileHeader().getSequence(mtSequence);
			if (samSequenceRecord == null)
				throw new RuntimeException("MT_SEQUENCE '" + mtSequence + "' is not found in sequence dictionary in " + INPUT.get(0).getAbsolutePath());
		}
	}

    /** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new SingleCellRnaSeqMetricsCollector().instanceMain(args));
	}
}
