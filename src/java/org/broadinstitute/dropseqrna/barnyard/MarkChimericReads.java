/*
 * MIT License
 *
 * Copyright 2021 Broad Institute
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

import com.google.common.base.CharMatcher;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.ChimericUmi.CHIMERIC_STRATEGY;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.OutputWriterUtil;
import org.broadinstitute.dropseqrna.utils.SamHeaderUtil;
import org.broadinstitute.dropseqrna.utils.StringInterner;
import org.broadinstitute.dropseqrna.utils.readiterators.GeneFunctionProcessor;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator;
import picard.cmdline.StandardOptionDefinitions;
import picard.nio.PicardHtsPath;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(summary = "Identify UMIs assigned to multiple genes in the same cell, then mark reads for such cell and UMI "
		+ "by adjusting the MAPQ, and/or generate a report for genes and UMIs in selected cells.", oneLineSummary = "Mark reads with UMIs that are assigned to multiple genes.", programGroup = DropSeq.class)
public class MarkChimericReads extends GeneFunctionCommandLineBase {
	private static final Log log = Log.getInstance(MarkChimericReads.class);
	
	private static final CHIMERIC_STRATEGY DEFAULT_STRATEGY = CHIMERIC_STRATEGY.RETAIN_MOST_SUPPORTED;
	public static final String CHIMERIC_COLUMN = "CHIMERIC";

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze. This argument can accept wildcards, or a file with the "
			+ "suffix .bam_list that contains the locations of multiple BAM files", minElements = 1)
	public List<PicardHtsPath> INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "SAM or BAM file with MAPQ of putative "
			+ "chimeric reads adjusted.", optional = true)
	public File OUTPUT;
	
	@Argument(doc="The strategy to use to detect and filter chimeric reads.  REMOVE_ALL marks all genes that share the same UMI on the same cell.  "
			+ "RETAIN_MOST_SUPPORTED retains the gene with the largest number of reads for a UMI/cell.  If RETAIN_MOST_SUPPORTED most suported gene is ambiguous,"
			+ "then all genes/UMI pairs are flagged as chimeric.", optional=false)
	public ChimericUmi.CHIMERIC_STRATEGY STRATEGY=DEFAULT_STRATEGY;

	@Argument(doc = "List of CELL_BARCODES to evaluate for chimeric reads.", optional=true)
	public File CELL_BC_FILE;

	@Argument(doc = "Tab-separated file with a row for each {CELL_BARCODE, MOLECULAR_BARCODE} found in selected cells, "
			+ "with the number of observations and chimeric status.  This file can end with .gz to be gzipped, which is recommended for larger data sets.", optional = true)
	public File OUTPUT_REPORT;

	@Argument(doc = "Produce report of number of unique {CELL_BARCODE, MOLECULAR_BARCODE}, and number marked.", shortName = StandardOptionDefinitions.METRICS_FILE_SHORT_NAME, optional = true)
	public File METRICS;

	@Argument(doc = "Mark as chimeric all reads in which the {CELL_BARCODE, MOLECULAR_BARCODE} is assigned to multiple genes.")
	public boolean MARK_UMI_REUSE = true;

	@Argument(doc = "Mark as chimeric reads with UMIs containing at least this many Ts.", optional = true)
	public Integer T_RICH_THRESHOLD;

	@Argument(doc = "Ignore reads with MAPQ < this when making chimeric determination.")
	public Integer READ_MQ = 10;

	@Argument(doc = "Tag in which to store original MAPQ, if read is marked as chimeric.  " + "Set to NULL to suppress tag creation.", optional = true)
	public String MAPQ_TAG = "YM";

	@Argument(doc = "Chimeric reads are tagged by setting MAPQ to this value.")
	public int CHIMERIC_MAPQ = 0;

	@Argument(doc = "The cell barcode tag.")
	public String CELL_BARCODE_TAG = "XC";

	@Argument(doc = "The molecular barcode tag.")
	public String MOLECULAR_BARCODE_TAG = "XM";

	@Override
	protected String[] customCommandLineValidation() {
		if (OUTPUT != null) {
			IOUtil.assertFileIsWritable(OUTPUT);
		}
		if (OUTPUT_REPORT != null) {
			IOUtil.assertFileIsWritable(OUTPUT_REPORT);
		}

		final ArrayList<String> list = new ArrayList<>(1);
		if (OUTPUT == null && OUTPUT_REPORT == null) {
			list.add("It does not make sense for neither OUTPUT nor OUTPUT_REPORT to be set.");
		}
		return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);
	}

	@Override
	protected int doWork() {
		INPUT = FileListParsingUtils.expandPicardHtsPathList(INPUT);

		// FIRST PASS: accumulate chimeric UMIs for each CBC, and optionally write report
		final Map<String, ChimericUmiCollection> chimerics = identifyChimericsAndWriteReport(MARK_UMI_REUSE, T_RICH_THRESHOLD);

		// PASS 2: Write BAM with chimeric reads marked.
		if (OUTPUT != null) {
			long numMarked = markReads(chimerics);
			log.info(numMarked + " reads marked chimeric.");
		}

		return 0;
	}

	private long markReads(final Map<String, ChimericUmiCollection> chimerics) {
		final GeneFunctionProcessor p = new GeneFunctionProcessor(GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG, false, STRAND_STRATEGY, LOCUS_FUNCTION_LIST, FUNCTIONAL_STRATEGY);
		
		long numMarked = 0;
		SamHeaderAndIterator headerAndIter = SamFileMergeUtil.mergeInputPaths(PicardHtsPath.toPaths(this.INPUT), false, SamReaderFactory.makeDefault());
		SamHeaderUtil.addPgRecord(headerAndIter.header, this);
		SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(headerAndIter.header, true, OUTPUT);
		ProgressLogger progLog = new ProgressLogger(log, 1000000, "marked");
		log.info("Marking chimeric UMIs in BAM");
		for (final SAMRecord rec : new IterableAdapter<>(headerAndIter.iterator)) {
			progLog.record(rec);
			// Parse the gene function tags, and get the interpreted tags for the given filters.
			// Using a deep copy of the record because otherwise processRead alters the record as a side effect.
			// TODO: should processRead  be making a deep copy internally to avoid the side effect?
			List<SAMRecord> processedRec =  p.processRead (rec.deepCopy());
			
			// only process the record for chimeric tags if the read can be interpreted with one consistent gene.
			if (processedRec.size()==1) {
				String geneName = processedRec.get(0) .getStringAttribute(GENE_NAME_TAG);
				String umi = rec.getStringAttribute(MOLECULAR_BARCODE_TAG);
				ChimericUmiCollection cuc = chimerics.get(rec.getStringAttribute(CELL_BARCODE_TAG));	
				if (cuc!=null && cuc.isChimericOrProblematic(umi, geneName)) {
					markChimeric(rec);
					++numMarked;
				}
			}						
			out.addAlignment(rec);
		}
		out.close();
		CloserUtil.close(headerAndIter.iterator);
		return numMarked;
	}

	private void markChimeric(final SAMRecord rec) {
		if (MAPQ_TAG != null) {
			rec.setAttribute(MAPQ_TAG, rec.getMappingQuality());
		}
		rec.setMappingQuality(CHIMERIC_MAPQ);
	}

	private static void writePerTranscriptHeader(final BufferedWriter out) {
		String[] header = { GatherMolecularBarcodeDistributionByGene.COLUMN_LABEL.CELL_BARCODE.toString(),
				GatherMolecularBarcodeDistributionByGene.COLUMN_LABEL.GENE.toString(),
				GatherMolecularBarcodeDistributionByGene.COLUMN_LABEL.MOLECULAR_BARCODE.toString(),
				GatherMolecularBarcodeDistributionByGene.COLUMN_LABEL.NUM_OBS.toString(), CHIMERIC_COLUMN};
		String h = StringUtils.join(header, "\t");
		OutputWriterUtil.writeResult(h, out);
	}
	
	private void writePerTranscriptStats(final String gene, final String cellBarcode, final ObjectCounter<String> counts, final ChimericUmiCollection chimerics,
			final BufferedWriter out, final MarkChimericReadMetrics metrics) {
				
		for (String key : counts.getKeys()) {
			int value = counts.getCountForKey(key);
			boolean chimeric = chimerics.isChimericOrProblematic(key, gene);
			String[] line = { cellBarcode, gene, key, Integer.toString(value), Boolean.toString(chimeric) };
			if (chimeric) {
				++metrics.NUM_MARKED_UMIS;
			}
			String h = StringUtils.join(line, "\t");
			OutputWriterUtil.writeResult(h, out);
		}
	}
	
	
	/**
     * @return Map with key=CBC, value=Set of chimeric UMIs for that CBC
     */
    private Map<String, ChimericUmiCollection> identifyChimericsAndWriteReport(boolean markUmiReuse, Integer tRichThreshold) {
        final BufferedWriter report_out;
        if (OUTPUT_REPORT != null) {
            IOUtil.assertFileIsWritable(OUTPUT_REPORT);
            report_out = IOUtil.openFileForBufferedWriting(OUTPUT_REPORT);

            writePerTranscriptHeader(report_out);
        } else {
            report_out = null;
        }
        
        // Set up the cell barcodes.  Can be null to try and repair the whole BAM.
        final Set<String> cellBarcodes=getCellBarcodes();       

		final List<Path> inputPaths = PicardHtsPath.toPaths(this.INPUT);
        PeekableIterator<UMICollection> umiIterator = new PeekableIterator<>(
                new UMIIterator.UMIIteratorBuilder(SamFileMergeUtil.mergeInputPaths(inputPaths, false),
                        GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG,
                        this.STRAND_STRATEGY, this.LOCUS_FUNCTION_LIST, this.FUNCTIONAL_STRATEGY, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG,
                        this.READ_MQ).setCellBarcodes(cellBarcodes).cellFirstSort(true).build());

        // Remember {CBC, UMI, Gene} pairs to be marked chimeric
        final Map<String, ChimericUmiCollection> chimerics = new HashMap<>();
        final MarkChimericReadMetrics metrics = new MarkChimericReadMetrics();

        // TODO: Is this premature optimization?  Genes will be reused across many cells so should be useful.
        StringInterner  geneStringCache = new StringInterner();
        
        // progress logger
        ProgressLogger progLog = new ProgressLogger(log, 10000, "cells chimeric marked");
        
        while (umiIterator.hasNext()) {
        	
            final String cellBarcode = umiIterator.peek().getCellBarcode();
            progLog.record(cellBarcode, 0);
            ChimericUmiCollection cuc = new ChimericUmiCollection(cellBarcode, this.STRATEGY);
            
            final List<UMICollection> umiCollectionsForCell = new ArrayList<>();
            
            // retain Alec's comments about slurping up data.
            while (umiIterator.hasNext() && cellBarcode.equals(umiIterator.peek().getCellBarcode())) {
            	UMICollection umiC= umiIterator.next();
            	// only add items if you want to find chimeras later.  If polyA only then this is false.
            	if (markUmiReuse) 
            		cuc.add(geneStringCache.intern(umiC.getGeneName()), umiC.getMolecularBarcodeCounts());
            		// cuc.add(umiC);
            	umiCollectionsForCell.add(umiC);
            }
            // finalize the collection.
            cuc.buildChimericUmis();
                        
            if (markUmiReuse) {                
                // Count number of UMIs for this CBC in case emitting metrics file.
                metrics.NUM_UMIS += cuc.getTotalUMIs();
                metrics.NUM_REUSED_UMIS += cuc.getTotalUMisChimeric();
            }
            
            if (tRichThreshold != null) {
                final Set<String> umisForCellBarcode = new HashSet<>();
                umiCollectionsForCell.forEach(umiCollection -> umisForCellBarcode.addAll(umiCollection.getMolecularBarcodeCounts().getKeys()));
                if (!markUmiReuse && METRICS != null) {
                    // Count number of UMIs for this CBC in case emitting metrics file, and not counted above
                    metrics.NUM_UMIS += umiCollectionsForCell.stream().mapToLong(umiCollection -> umiCollection.getMolecularBarcodeCounts().getSize()).sum();
                }
                final CharMatcher tMatcher = CharMatcher.is('T');
                final Collection<String> polyTUmis = umisForCellBarcode.stream().filter(umi -> tMatcher.countIn(umi) >= tRichThreshold).collect(Collectors.toList());
                metrics.NUM_T_RICH_UMIS += polyTUmis.size();
                polyTUmis.forEach(cuc::registerProblematicUmi);
            }
                        
            chimerics.put(cellBarcode, cuc);
            if (report_out != null) {
                for (final UMICollection batch : umiCollectionsForCell) {
                	
                    writePerTranscriptStats(batch.getGeneName(), batch.getCellBarcode(),
                            batch.getMolecularBarcodeCounts(), cuc, report_out, metrics);
                }
            }
        }
        CloserUtil.close(umiIterator);
        
        log.info("Total cell/UMIs observed [" +metrics.NUM_UMIS+"] Marked Chimeric [" + metrics.NUM_MARKED_UMIS+"] % ["+ String.format("%.2f%%",metrics.getPercentMarked())+"]");
        
        try {
            if (report_out != null) {
                report_out.close();
            }
        } catch (IOException e) {
            throw new RuntimeIOException("Exception writing " + OUTPUT_REPORT, e);
        }

        if (METRICS != null) {
            final MetricsFile<MarkChimericReadMetrics, Integer> metricsFile = new MetricsFile<>();
            metricsFile.addMetric(metrics);
            metricsFile.write(METRICS);
        }

        return chimerics;
    }
    
    private Set<String> getCellBarcodes () {
		if (CELL_BC_FILE==null)  {
			log.info("No cell barcode file - will process all cell barcodes");
			return null;
        }        	        
         	
		Set<String> cellBarcodes=new HashSet<>(ParseBarcodeFile.readCellBarcodeFile(CELL_BC_FILE));
    	log.info("Found " + cellBarcodes.size()+ " cell barcodes in file");
    	return cellBarcodes;
	}
	
	public static class MarkChimericReadMetrics extends MetricBase {
		public long NUM_UMIS;
		public long NUM_MARKED_UMIS;
		public long NUM_T_RICH_UMIS;
		public long NUM_REUSED_UMIS;
		
		public double getPercentMarked () {
			return ((double) NUM_MARKED_UMIS / (double) NUM_UMIS)*100;
		}

		/**
		 * Adds the non-calculated metrics (which is all of them)
		 */
		public void merge(final MarkChimericReadMetrics metric) {
			this.NUM_UMIS += metric.NUM_UMIS;
			this.NUM_MARKED_UMIS += metric.NUM_MARKED_UMIS;
			this.NUM_T_RICH_UMIS += metric.NUM_T_RICH_UMIS;
			this.NUM_REUSED_UMIS += metric.NUM_REUSED_UMIS;
		}
	}

	/** Stock main method. */
	public static void main(final String[] args) {
		new MarkChimericReads().instanceMainWithExit(args);
	}
}
