package org.broadinstitute.dropseqrna.metrics;

import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.GeneFunctionCommandLineBase;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.*;
import org.broadinstitute.dropseqrna.utils.readiterators.*;
import org.broadinstitute.dropseqrna.utils.readiterators.umisort.ReadIntervalSummary;
import picard.cmdline.StandardOptionDefinitions;
import picard.nio.PicardHtsPath;

import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

@CommandLineProgramProperties(
        summary = "For a cell barcode/gene/molecular barcode, finds the min/max positions of the reads in genomic coordinates.  The report contains the following columns:\n" +
				"Cell_Barcode, Gene, Molecular_Barcode: there is a single row for each unique combination of these.\n" +
				"Num_Reads: total reads for this {cell barcode, gene, molecular barcode}\n" +
				"Num_Gapped_Reads: number of reads for this row that have at least one gap\n" +
				"Contig: contig containing the gene\n" +
				"Position_Min: smallest 1-based alignment start for all reads for this row\n" +
				"Position_Max: largest 1-based inclusive alignment end for all reads for this row\n" +
				"Read_Positive_Strand: true if reads are on positive strand\n" +
				"Gaps: comma-separated list of gaps, where each gap is <1-based position before the gap start>:<1-based position after the gap end>:<number of reads with this gap>",
        oneLineSummary = "Find read intervals for each UMI",
        programGroup = DropSeq.class
)
public class GatherUMIReadIntervals extends GeneFunctionCommandLineBase {
	
	private static final Log log = Log.getInstance(GatherUMIReadIntervals.class);
	
	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze. This argument can accept wildcards, or a file with the suffix .bam_list that contains the locations of multiple BAM files", minElements = 1)
	public List<PicardHtsPath> INPUT;
	
	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output file of read intervals per UMI. This supports zipped formats like gz and bz2.")
	public File OUTPUT;

	@Argument(doc="An interval string, defined as contig:start-end.  If supplied only reads that overlap this interval are reported.", optional = true)
	public String INTERVAL;

	@Argument(doc="An interval file, which contains a list of intervals to query.  If supplied only reads that overlap these intervals are reported. " +
			"The interval file is in Picard-style format - see https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists for details", optional = true)
	public File INTERVAL_LIST;

	@Argument(doc = "The cell barcode tag.  If there are no reads with this tag, the program will assume that all reads belong to the same cell and process in single sample mode.")
	public String CELL_BARCODE_TAG = "XC";

	@Argument(doc = "The molecular barcode tag.")
	public String MOLECULAR_BARCODE_TAG = "XM";

	@Argument(doc = "The map quality of the read to be included.")
	public Integer READ_MQ = 10;

	@Argument(doc = "Filters data to a subset of cell barcodes in the input BAM.  If left unset, all cell barcodes are analyzed.  The file has 1 column with no header.", optional=true)
	public File CELL_BC_FILE = null;

	@Argument(doc="The edit distance that molecular barcodes should be combined at within a gene.")
	public int EDIT_DISTANCE=0;

	@Argument(doc = "if true, records tagged with multiple genes be double counted, once for each gene.  " +
			"If false, reads that map to multiple genes are ignored.")
	public boolean ASSIGN_READS_TO_ALL_GENES=true;

	@Argument(doc = "Number of threads used to read and decode the input BAM.  " +
			"A single BAM speeds up as threads are added, up to about 2 threads.")
	public int NUM_THREADS = 1;

	@Override
	protected int doWork() {
		this.INPUT = FileListParsingUtils.expandPicardHtsPathList(INPUT);
		IOUtil.assertFileIsWritable(this.OUTPUT);
		if (CELL_BC_FILE!=null) IOUtil.assertFileIsReadable(this.CELL_BC_FILE);
	
		BufferedWriter out = IOUtil.openFileForBufferedWriting(OUTPUT);
		writePerUMIStatsHeader(out);
		
		List<String> cellBarcodes=getCellBarcodes(this.CELL_BC_FILE);
		final SamHeaderAndIterator headerAndIterator =
				SamFileMergeUtil.mergeInputPaths(PicardHtsPath.toPaths(this.INPUT), false,
						SamReaderFactory.makeDefault(), this.NUM_THREADS);

		IntervalList intervals = getInputIntervals(headerAndIterator.header);

		// Cell first, so the rows come out ordered by cell, gene, then UMI.
		final UMIIterator umiIterator = new UMIIterator.UMIIteratorBuilder(headerAndIterator, this.GENE_NAME_TAG,
				this.GENE_STRAND_TAG, this.GENE_FUNCTION_TAG, this.STRAND_STRATEGY, this.LOCUS_FUNCTION_LIST,
				this.FUNCTIONAL_STRATEGY, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, this.READ_MQ)
				.assignReadsToAllGenes(ASSIGN_READS_TO_ALL_GENES).setCellBarcodes(cellBarcodes).setIntervals(intervals)
				.cellFirstSort(true).collectReadIntervals(true).build();

		int counter=0;
		for (final UMICollection c : new IterableAdapter<>(umiIterator)) {
			if (EDIT_DISTANCE > 0)
				c.collapseThisByEditDistance(EDIT_DISTANCE);
			final List<String> umis = new ArrayList<>(c.getMolecularBarcodes());
			Collections.sort(umis);
			for (final String umi : umis) {
				if (counter%1000000==0)
					log.info("Processed [" + Math.round(counter/1e6) + "] Million UMIs");
				final UmiReadInterval i = new UmiReadInterval(c.getCellBarcode(), c.getGeneName(), umi,
						c.getMolecularBarcodeCounts().getCountForKey(umi), c.getReadInterval(umi));
				writePerUMIStats(i, out);
				counter++;
			}
		}
		log.info("Finished Scan. Processed [" + counter + "] UMIs");

		umiIterator.close();
		CloserUtil.close(out);
		return 0;

	}

	@Override
	protected String[] customCommandLineValidation() {
		final ArrayList<String> list = new ArrayList<>();
		if (this.EDIT_DISTANCE < 0) {
			list.add("EDIT_DISTANCE cannot be negative.");
		}
		if (this.NUM_THREADS < 1) {
			list.add("NUM_THREADS must be >= 1");
		}
		return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);
	}

	/**
	 * If a string interval or interval file is provided, parse and generate an IntervalList.
	 * Otherwise emit null.
	 * @param header The BAM header.  The interval contig should match up to the BAM contig.
	 * @return An interval list or null.
	 */
	private IntervalList getInputIntervals (SAMFileHeader header) {
		if (this.INTERVAL!=null) {
			Interval i= parseInterval(this.INTERVAL);
			IntervalList result = new IntervalList(header);
			result.add(i);
			return result;
		}
		if (this.INTERVAL_LIST!=null) {
			IOUtil.assertFileIsReadable(this.INTERVAL_LIST);
			return IntervalList.fromFile(this.INTERVAL_LIST);
		}
		return null;
	}

	private Interval parseInterval(String intervalString) {
		String[] parts = intervalString.split(":");
		if (parts.length != 2) {
			throw new IllegalArgumentException("Invalid interval format: " + intervalString);
		}

		String contig = parts[0];
		String[] range = parts[1].split("-");
		if (range.length != 2) {
			throw new IllegalArgumentException("Invalid interval format: " + intervalString);
		}

		int start = Integer.parseInt(range[0]);
		int end = Integer.parseInt(range[1]);

		return new Interval(contig, start, end);
	}
	
	static void writePerUMIStats (UmiReadInterval interval, final BufferedWriter out) {
		OutputWriterUtil.writeResult(interval.toTsv(), out);

	}

	static void writePerUMIStatsHeader(final BufferedWriter out) {
		String h = StringUtils.join(UmiReadInterval.COLUMN_HEADERS, "\t");
		OutputWriterUtil.writeResult(h, out);
	}
	
	private List<String> getCellBarcodes (File cellBCFile) {
		if (cellBCFile!=null) {
			List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(cellBCFile);
			log.info("Found " + cellBarcodes.size()+ " cell barcodes in file");
			return (cellBarcodes);
		}
		return null;

	}
	
	
	static class UmiReadInterval {
		public static final String CELL_BARCODE_HEADER = "Cell_Barcode";
		public static final String GENE_HEADER = "Gene";
		public static final String MOLECULAR_BARCODE_HEADER = "Molecular_Barcode";
		public static final String NUM_READS_HEADER = "Num_Reads";
		public static final String NUM_GAPPED_READS_HEADER = "Num_Gapped_Reads";
		public static final String CONTIG_HEADER = "Contig";
		public static final String POSITION_MIN_HEADER = "Position_Min";
		public static final String POSITION_MAX_HEADER = "Position_Max";
		public static final String READ_POSITIVE_STRAND_HEADER = "Read_Positive_Strand";
		public static final String GAPS_HEADER = "Gaps";
		static String [] COLUMN_HEADERS = {CELL_BARCODE_HEADER,
				GENE_HEADER, MOLECULAR_BARCODE_HEADER, NUM_READS_HEADER, NUM_GAPPED_READS_HEADER, CONTIG_HEADER,
				POSITION_MIN_HEADER, POSITION_MAX_HEADER, READ_POSITIVE_STRAND_HEADER, GAPS_HEADER};

		final String CELL_BARCODE;
		final String GENE_NAME;
		final String MOLECULAR_BARCODE;
		final int NUM_READS;
		final ReadIntervalSummary INTERVAL;

		UmiReadInterval(final String cellBarcode, final String geneName, final String molecularBarcode,
						final int numReads, final ReadIntervalSummary interval) {
			this.CELL_BARCODE = cellBarcode;
			this.GENE_NAME = geneName;
			this.MOLECULAR_BARCODE = molecularBarcode;
			this.NUM_READS = numReads;
			this.INTERVAL = interval;
		}

		@Override
		public String toString() {
			return "UmiReadInterval [CELL_BARCODE=" + CELL_BARCODE + ", GENE_NAME=" + GENE_NAME + ", MOLECULAR_BARCODE="
					+ MOLECULAR_BARCODE + ", NUM_READS=" + NUM_READS + ", INTERVAL=" + INTERVAL + "]";
		}

		/** @return each gap as start:end:count, ordered by position, comma separated. */
		private String makeGapsString() {
			final List<String> gapStrings = new ArrayList<>(INTERVAL.getGaps().size());
			for (Map.Entry<ReadIntervalSummary.Gap, Integer> entry : INTERVAL.getGaps().entrySet()) {
				final ReadIntervalSummary.Gap gap = entry.getKey();
				gapStrings.add(gap.start() + ":" + gap.end() + ":" + entry.getValue());
			}
			return StringUtils.join(gapStrings, ",");
		}

		public String toTsv() {
			String [] body = {CELL_BARCODE, GENE_NAME, MOLECULAR_BARCODE, Integer.toString(NUM_READS),
					Integer.toString(INTERVAL.getNumGappedReads()), INTERVAL.getContig(),
					Integer.toString(INTERVAL.getPositionMin()), Integer.toString(INTERVAL.getPositionMax()),
					Boolean.toString(INTERVAL.isPositiveStrand()), makeGapsString()};
			return StringUtils.join(body, "\t");
		}
	}
	
	
	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new GatherUMIReadIntervals().instanceMain(args));
	}
	
}
