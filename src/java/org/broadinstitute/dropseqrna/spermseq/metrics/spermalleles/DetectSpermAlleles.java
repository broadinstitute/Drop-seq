package org.broadinstitute.dropseqrna.spermseq.metrics.spermalleles;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.GeneFunctionCommandLineBase;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.DigitalAlleleCounts;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.DigitalAlleleCountsIterator;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.MultiCellDigitalAlleleCounts;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.MultiCellDigitalAlleleCountsIterator;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPUMIBasePileupIterator;
import org.broadinstitute.dropseqrna.cmdline.SpermSeq;
import org.broadinstitute.dropseqrna.priv.barnyard.digitalallelecounts.GatherDigitalAlleleCounts;
import org.broadinstitute.dropseqrna.priv.barnyard.digitalallelecounts.sampleassignment.SortOrder;
import org.broadinstitute.dropseqrna.utils.CountChangingIteratorWrapper;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import picard.annotation.LocusFunction;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;


@CommandLineProgramProperties(
        summary = "Detect which SNPs are present in each sperm cell (SNPs typed are from the donor’s genome, provided via interval file). For each SNP present, report the number of reads and UMIs containing each base (A, C, T, G, or N).",
        oneLineSummary = "Detect which alleles of which SNPs are present in each sperm cell”",
        programGroup = SpermSeq.class
)
public class DetectSpermAlleles extends CommandLineProgram {
	
		private static final Log log = Log.getInstance(DetectSpermAlleles.class);

		@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.", optional=false)
		public File INPUT;

		@Argument(doc = "The list of snp intervals to gather SpermSeq allele counts on.  This file is in Interval format - tab seperated with fields: chr start end strand name", optional=false)
		public File INTERVALS;

		@Argument(doc="File contain cell barcodes to process. The file has 1 column with no header.", optional=false)
		public File CELL_BC_FILE=null;

		@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file of SpermSeq alleles.  One SNP/Sperm cell per row. This supports zipped formats like gz and bz2.")
		public File OUTPUT;
		
		@Argument(doc="The cell barcode tag.  If there are no reads with this tag, the program will assume that all reads belong to the same cell and process in single sample mode.")
		public String CELL_BARCODE_TAG="XC";

		@Argument(doc="The molecular barcode tag.")
		public String MOLECULAR_BARCODE_TAG="XM";

		@Argument(doc="The edit distance that molecular barcodes should be combined at within a SNP.")
		public Integer EDIT_DISTANCE=1;

		@Argument(doc="The map quality of the read to be included.")
		public Integer READ_MQ=10;

		@Argument(doc="The base quality of the read bases to be included.")
		public Integer BASE_QUALITY=10;
				
		// @Argument (doc="Remove UMIs below this purity threshold.  A UMI's purity is determined as the number of reads of the most common base divided by the total number of reads.")
		private Double UMI_PURITY_THRESHOLD=1.0;

		private final boolean AUTO_FLUSH_OUTPUTS=true;
		private final String SNP_TAG = "YS";
		private final char [] BASES = {'A', 'C', 'G', 'T', 'N'};
		private final int PROGRESS_RATE=100000;

		@Override
		protected int doWork() {
			// validation
			IOUtil.assertFileIsReadable(INPUT);
			IOUtil.assertFileIsWritable(OUTPUT);
			IOUtil.assertFileIsReadable(this.CELL_BC_FILE);
			
			List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(this.CELL_BC_FILE);
			log.info("Using " + cellBarcodes.size() + " cells in analysis");
			IntervalList snpIntervals = IntervalList.fromFile(INTERVALS);
			log.info("Using " + snpIntervals.getIntervals().size() + " SNP intervals in analysis");
			PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(OUTPUT));
			writeHeader(out);

			SamReader reader = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE).open(this.INPUT);			
			
			SNPUMIBasePileupIterator sbpi = getIter(reader, snpIntervals, cellBarcodes);
			
			MultiCellDigitalAlleleCountsIterator multiIter = new MultiCellDigitalAlleleCountsIterator(new DigitalAlleleCountsIterator(sbpi, BASE_QUALITY));

			// sort cell barcodes alphabetically for output.
			Collections.sort(cellBarcodes);
			@SuppressWarnings("unused")
			int counter=0;
			while (multiIter.hasNext()) {
				MultiCellDigitalAlleleCounts mcdac = multiIter.next();
				processMCDAC(cellBarcodes, mcdac, out, AUTO_FLUSH_OUTPUTS);
				counter++;
				if (counter%PROGRESS_RATE==0) log.info("Processed " + counter + " SNPs");
			}
			log.info("Processed " + counter +" total SNPs");
			out.close();
			multiIter.close();
			return 0;
		}

		/**
		 * Takes all the cells for a SNP/GENE and writes out their summary statistics.
		 * @param cellBarcodes The list of cells to try and write out - this provides output ordering once you're past the SNP/GENE level.
		 * @param mcdac
		 * @param out The per-cell output stream
		 * @param outSum The meta analysis output stream
		 */
		private void processMCDAC (final List<String> cellBarcodes, final MultiCellDigitalAlleleCounts mcdac, final PrintStream out, final boolean autoFlush) {

			// collapses all DAC UMIs.
			mcdac.collapseDACs(this.EDIT_DISTANCE);
			mcdac.filterDataByUMIPurity(UMI_PURITY_THRESHOLD);
			// write meta analysis stats for non-empty objects.
			DigitalAlleleCounts metaDAC = mcdac.getMetaAnalysis();
			// if the meta analysis has no information (all reads filtered out?)  return.
			if (metaDAC.isEmpty())
				return;
								
			// loop over cells and write stats for non-empty objects.
			for (String cell: cellBarcodes) {
				DigitalAlleleCounts dac = mcdac.getDigitalAlleleCounts(cell);
				if (dac!=null && !dac.isEmpty()) {					
					writeStats(dac, out);
					if (autoFlush) out.flush();
				}
			}
		}

		private void writeStats (final DigitalAlleleCounts dac, final PrintStream out) {
			List<String> line = new ArrayList<>();
			line.add(dac.getSnpInterval().getContig());
			line.add(Integer.toString(dac.getSnpInterval().getStart()));			
			line.add(dac.getCell());			
			ObjectCounter<Character> readCounts = dac.getReadCounts();
			for (Character base: BASES)
				line.add(Integer.toString(readCounts.getCountForKey(base)));

			ObjectCounter<Character> umiCounts = dac.getUMIAlleleCount();
			for (Character base: BASES)
				line.add(Integer.toString(umiCounts.getCountForKey(base)));
			
			String h = StringUtils.join(line, "\t");
			out.println(h);
		}

		private void writeHeader(final PrintStream out) {
			List<String> header = new ArrayList<>();

			header.add("chr");
			header.add("pos");
			header.add("cell");
			for (char c : BASES)
				header.add("r"+c);
			for (char c : BASES)
				header.add(c+"");
			
			String h = StringUtils.join(header, "\t");
			out.println(h);
		}

		private SNPUMIBasePileupIterator getIter (SamReader reader, IntervalList snpIntervals, List <String> cellBarcodes) {
			
			Iterator<SAMRecord> iter = reader.iterator();
			
			// pre-process iterator to add the "gene" tag and the gene strand.
			TagAddingIterator taggedIter = new TagAddingIterator(iter);						
			SamHeaderAndIterator headerAndIter = new SamHeaderAndIterator(reader.getFileHeader(), taggedIter);			
			
			// some hard coded params.	
			
			SNPUMIBasePileupIterator sbpi = new SNPUMIBasePileupIterator(
					headerAndIter, snpIntervals, GeneFunctionCommandLineBase.DEFAULT_GENE_NAME_TAG, GeneFunctionCommandLineBase.DEFAULT_GENE_STRAND_TAG, GeneFunctionCommandLineBase.DEFAULT_GENE_FUNCTION_TAG,
					GeneFunctionCommandLineBase.DEFAULT_LOCUS_FUNCTION_LIST, StrandStrategy.BOTH, this.CELL_BARCODE_TAG,
					this.MOLECULAR_BARCODE_TAG, this.SNP_TAG, null, this.READ_MQ, true, cellBarcodes, SortOrder.SNP_GENE);
			
			return sbpi;
			
		}
		
		/**
		 * Adds the gene function annotations to reads on the fly.
		 * @author nemesh
		 *
		 */
		private class TagAddingIterator extends CountChangingIteratorWrapper<SAMRecord>  {
			
			protected TagAddingIterator(Iterator<SAMRecord> underlyingIterator) {
				super(underlyingIterator);
			}

			@Override
			protected void processRecord(SAMRecord rec) {
				rec.setAttribute(GeneFunctionCommandLineBase.DEFAULT_GENE_NAME_TAG, rec.getContig());
				rec.setAttribute(GeneFunctionCommandLineBase.DEFAULT_GENE_STRAND_TAG, "-");
				String coding = LocusFunction.CODING.toString();
				rec.setAttribute(GeneFunctionCommandLineBase.DEFAULT_GENE_FUNCTION_TAG, LocusFunction.CODING.toString());
				queueRecordForOutput(rec);					
			}
						
		}
		
		/** Stock main method. */
		public static void main(final String[] args) {
			System.exit(new DetectSpermAlleles().instanceMain(args));
		}
}


