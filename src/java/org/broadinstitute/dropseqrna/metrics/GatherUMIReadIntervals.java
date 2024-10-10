package org.broadinstitute.dropseqrna.metrics;

import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

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

import picard.annotation.LocusFunction;
import picard.cmdline.StandardOptionDefinitions;
import picard.nio.PicardHtsPath;
import picard.util.TabbedTextFileWithHeaderParser;

@CommandLineProgramProperties(
        summary = "For a cell/gene/umi, finds the min/max positions of the reads in genomic coordinates.",
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


	@Override
	protected int doWork() {
		this.INPUT = FileListParsingUtils.expandPicardHtsPathList(INPUT);
		IOUtil.assertFileIsWritable(this.OUTPUT);
		if (CELL_BC_FILE!=null) IOUtil.assertFileIsReadable(this.CELL_BC_FILE);
	
		BufferedWriter out = IOUtil.openFileForBufferedWriting(OUTPUT);
		writePerUMIStatsHeader(out);
		
		List<String> cellBarcodes=getCellBarcodes(this.CELL_BC_FILE);
		final SamHeaderAndIterator headerAndIterator =
				SamFileMergeUtil.mergeInputPaths(PicardHtsPath.toPaths(this.INPUT), false);

		IntervalList intervals = getInputIntervals(headerAndIterator.header);

		GroupingIterator<SAMRecord> iter = prepareIter(headerAndIterator, this.CELL_BARCODE_TAG, this.GENE_NAME_TAG, this.MOLECULAR_BARCODE_TAG, 
				this.GENE_STRAND_TAG, this.GENE_FUNCTION_TAG, this.STRAND_STRATEGY, this.LOCUS_FUNCTION_LIST, this.READ_MQ, cellBarcodes, intervals);
		
		int counter=0;
		while (iter.hasNext()) {			
			if (counter%1000000==0)
				log.info("Processed [" + Math.round(counter/1e6) + "] Million UMIs");
			List<SAMRecord> recs = iter.next();
			UmiReadInterval i = getInterval(recs, this.CELL_BARCODE_TAG, this.GENE_NAME_TAG, this.MOLECULAR_BARCODE_TAG);
			writePerUMIStats(i, out);
			counter++;
		}
		log.info("Finished Scan. Processed [" + counter + "] UMIs");
		
		headerAndIterator.iterator.close();
		CloserUtil.close(out);		
		return 0;
		
	}

	@Override
	protected String[] customCommandLineValidation() {
		final ArrayList<String> list = new ArrayList<>();
		if (this.EDIT_DISTANCE < 0) {
			list.add("EDIT_DISTANCE cannot be negative.");
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
	
	private UmiReadInterval getInterval (List<SAMRecord> recs, String cellBarcodeTag, String geneTag, String molecularBarcodeTag) {
		UmiReadInterval result = new UmiReadInterval();
		result.CELL_BARCODE=recs.getFirst().getStringAttribute(cellBarcodeTag);
		result.GENE_NAME=recs.getFirst().getStringAttribute(geneTag);
		result.MOLECULAR_BARCODE=recs.getFirst().getStringAttribute(molecularBarcodeTag);
		result.NUM_READS=recs.size();
		result.POSITIVE_STRAND=!recs.getFirst().getReadNegativeStrandFlag();
		result.CONTIG=recs.getFirst().getContig();
		
		// gather up the alignment positions for the min/max
		List<Integer> posList = new ArrayList<>();		
		
		for (SAMRecord r: recs) {
			// get the min/max, count of reads with cigar string N
			if (hasGap(r)) result.NUM_SPLICED_READS++;	
			posList.add(r.getAlignmentStart());
			posList.add(r.getAlignmentEnd());
		}
		
		result.POSITION_MIN=posList.stream().mapToInt(x->x).min().getAsInt();
		result.POSITION_MAX=posList.stream().mapToInt(x->x).max().getAsInt();
				
		return (result);
	}
	
	private boolean hasGap (final SAMRecord r) {
		Cigar c = r.getCigar();
		for (CigarElement ce: c.getCigarElements())
			if (ce.getOperator()==CigarOperator.N)
				return true;
		return false;
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
	
	
	private GroupingIterator<SAMRecord> prepareIter (SamHeaderAndIterator headerAndIterator, String cellBarcodeTag, String geneTag, String molecularBarcodeTag,
													 String geneStrandTag, String geneFunctionTag, StrandStrategy strandStrategy,  List<LocusFunction> acceptedLociFunctions,
													 Integer readMQ, Collection <String> cellBarcodes, IntervalList intervals) {

		final StringTagComparator cellBarcodeTagComparator = new StringTagComparator(cellBarcodeTag);
        final StringTagComparator geneExonTagComparator = new StringTagComparator(geneTag);
        final StringTagComparator umiComparator = new StringTagComparator (molecularBarcodeTag);
        
        final MultiComparator<SAMRecord> multiComparator = new MultiComparator<>(cellBarcodeTagComparator, geneExonTagComparator, umiComparator);

		final UMICollectionIteratorAdapter umiIterator = new UMICollectionIteratorAdapter(
				new UMIIterator.UMIIteratorBuilder(headerAndIterator, geneTag, geneStrandTag, geneFunctionTag,
				strandStrategy, acceptedLociFunctions, this.FUNCTIONAL_STRATEGY, cellBarcodeTag, molecularBarcodeTag,
				readMQ).assignReadsToAllGenes(ASSIGN_READS_TO_ALL_GENES).setCellBarcodes(cellBarcodes).retainReads(true).
						setIntervals(intervals).build());

		// Filter/assign reads based on functional annotations
		GeneFunctionIteratorWrapper gfteratorWrapper = new GeneFunctionIteratorWrapper(umiIterator, geneTag,
				geneStrandTag, geneFunctionTag, true, strandStrategy, acceptedLociFunctions, this.FUNCTIONAL_STRATEGY);


        CloseableIterator<SAMRecord> sortedAlignmentIterator = SamRecordSortingIteratorFactory.create(
                headerAndIterator.header, gfteratorWrapper, multiComparator, null);

        GroupingIterator<SAMRecord> finalIter = new GroupingIterator<>(sortedAlignmentIterator, multiComparator);
        return (finalIter);

	}

	/**
	 * Convert a stream of UMICollections into a stream of SAMRecords, with the reads for each UMI and sorted
	 * by UMI
	 */
	class UMICollectionIteratorAdapter
			implements CloseableIterator<SAMRecord> {
		private final UMIIterator umiCollectionIterator;
		private Iterator<String> umiIterator = Collections.emptyIterator();
		private Iterator<SAMRecord> readIterator = Collections.emptyIterator();
		private UMICollection currentCollection = null;
		public UMICollectionIteratorAdapter(final UMIIterator umiCollectionIterator) {
			this.umiCollectionIterator = umiCollectionIterator;
		}

		@Override
		public void close() {
			this.umiCollectionIterator.close();
		}

		@Override
		public boolean hasNext() {
			return readIterator.hasNext() || umiIterator.hasNext() || umiCollectionIterator.hasNext();
		}

		@Override
		public SAMRecord next() {
			if (readIterator.hasNext()) return readIterator.next();
			if (umiIterator.hasNext()) {
				String umi = umiIterator.next();
				readIterator = this.currentCollection.getReads(umi).iterator();
				if (!readIterator.hasNext()) throw new IllegalStateException("No reads found for UMI [" + umi + "]");
				return next();
			}
			if (!umiCollectionIterator.hasNext()) throw new IllegalStateException("Should not call next() if !hasNext()");
			currentCollection = umiCollectionIterator.next();
			if (EDIT_DISTANCE > 0) {
				currentCollection.collapseThisByEditDistance(EDIT_DISTANCE, MOLECULAR_BARCODE_TAG);
			}
			final List<String> umis = new ArrayList<>(currentCollection.getMolecularBarcodes());
			Collections.sort(umis);
			umiIterator = umis.iterator();
			if (!umiIterator.hasNext()) throw new IllegalStateException("No UMIs found for collection");
			return next();
		}
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
		static String [] COLUMN_HEADERS = {CELL_BARCODE_HEADER,
				GENE_HEADER, MOLECULAR_BARCODE_HEADER, NUM_READS_HEADER, NUM_GAPPED_READS_HEADER, CONTIG_HEADER,
				POSITION_MIN_HEADER, POSITION_MAX_HEADER, READ_POSITIVE_STRAND_HEADER};

		String CELL_BARCODE;
		String GENE_NAME;
		String MOLECULAR_BARCODE;
		int NUM_READS;
		int NUM_SPLICED_READS=0;
		String CONTIG;
		int POSITION_MIN;
		int POSITION_MAX;
		boolean POSITIVE_STRAND;

		public UmiReadInterval() {
		}

		public UmiReadInterval(final TabbedTextFileWithHeaderParser.Row row) {
			this.CELL_BARCODE = row.getField(CELL_BARCODE_HEADER);
			this.GENE_NAME = row.getField(GENE_HEADER);
			this.MOLECULAR_BARCODE = row.getField(MOLECULAR_BARCODE_HEADER);
			this.NUM_READS = row.getIntegerField(NUM_READS_HEADER);
			this.NUM_SPLICED_READS = row.getIntegerField(NUM_GAPPED_READS_HEADER);
			this.CONTIG = row.getField(CONTIG_HEADER);
			this.POSITION_MIN = row.getIntegerField(POSITION_MIN_HEADER);
			this.POSITION_MAX = row.getIntegerField(POSITION_MAX_HEADER);
			this.POSITIVE_STRAND = Boolean.parseBoolean(row.getField(READ_POSITIVE_STRAND_HEADER));
		}

		@Override
		public String toString() {
			return "UmiReadInterval [CELL_BARCODE=" + CELL_BARCODE + ", GENE_NAME=" + GENE_NAME + ", MOLECULAR_BARCODE="
					+ MOLECULAR_BARCODE + ", NUM_READS=" + NUM_READS + ", NUM_SPLICED_READS=" + NUM_SPLICED_READS
					+ ", CONTIG=" + CONTIG + ", POSITION_MIN=" + POSITION_MIN + ", POSITION_MAX=" + POSITION_MAX
					+ ", POSITIVE_STRAND=" + POSITIVE_STRAND + "]";
		}		

		public String toTsv() {
			String [] body = {CELL_BARCODE, GENE_NAME, MOLECULAR_BARCODE, Integer.toString(NUM_READS), Integer.toString(NUM_SPLICED_READS),
					CONTIG, Integer.toString(POSITION_MIN), Integer.toString(POSITION_MAX), Boolean.toString(POSITIVE_STRAND)};
			return StringUtils.join(body, "\t");
		}
	}
	
	
	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new GatherUMIReadIntervals().instanceMain(args));
	}
	
}
