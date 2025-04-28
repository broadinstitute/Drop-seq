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
import picard.annotation.LocusFunction;
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
			final List<AlignmentGap> gaps = getGaps(r);
			if (!gaps.isEmpty()) {
				result.NUM_SPLICED_READS++;
				for (final AlignmentGap gap: gaps) {
					if (result.GAPS.containsKey(gap)) {
						result.GAPS.put(gap, result.GAPS.get(gap) + 1);
					} else {
						result.GAPS.put(gap, 1);
					}
				}
			}
			posList.add(r.getAlignmentStart());
			posList.add(r.getAlignmentEnd());
		}
		
		result.POSITION_MIN=posList.stream().mapToInt(x->x).min().getAsInt();
		result.POSITION_MAX=posList.stream().mapToInt(x->x).max().getAsInt();
				
		return (result);
	}
	
	private List<Integer> getIndicesOfAlignmentBlocksPrecedingGaps(final SAMRecord r) {
		ArrayList<Integer> ret = null;
		List<CigarElement> cigarElements = r.getCigar().getCigarElements();
		int alignmentBlockIndex = -1;
		for (final CigarElement ce: cigarElements) {
			if (ce.getOperator() == CigarOperator.MATCH_OR_MISMATCH) {
				++alignmentBlockIndex;
			} else if (ce.getOperator() == CigarOperator.SKIPPED_REGION) {
				if (ret == null) {
					ret = new ArrayList<>();
				}
				ret.add(alignmentBlockIndex);
			}
		}
		if (ret == null) {
			return Collections.emptyList();
		}
		return ret;
	}

	private List<AlignmentGap> getGaps (final SAMRecord r) {
		final List<Integer> alignmentBlocksPrecedingGaps = getIndicesOfAlignmentBlocksPrecedingGaps(r);
		if (alignmentBlocksPrecedingGaps.isEmpty()) {
			return Collections.emptyList();
		} else {
			List<AlignmentGap> ret = new ArrayList<>(alignmentBlocksPrecedingGaps.size());
			final List<AlignmentBlock> alignmentBlocks = r.getAlignmentBlocks();
			for (int alignmentBlockIndex : alignmentBlocksPrecedingGaps) {
				final AlignmentBlock precedingAlignmentBlock = alignmentBlocks.get(alignmentBlockIndex);
				if (alignmentBlockIndex == alignmentBlocks.size() - 1) {
					throw new IllegalStateException("Strange CIGAR: " + r.getCigarString() +
							" for read " + r.getReadName());
				}
				ret.add(new AlignmentGap(precedingAlignmentBlock.getReferenceStart() + precedingAlignmentBlock.getLength() - 1,
						alignmentBlocks.get(alignmentBlockIndex + 1).getReferenceStart()));
			}
			return ret;
		}
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

	record AlignmentGap(int start, int end) {

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

		String CELL_BARCODE;
		String GENE_NAME;
		String MOLECULAR_BARCODE;
		int NUM_READS;
		int NUM_SPLICED_READS=0;
		String CONTIG;
		int POSITION_MIN;
		int POSITION_MAX;
		boolean POSITIVE_STRAND;
		Map<AlignmentGap, Integer> GAPS = new HashMap<>();

		public UmiReadInterval() {
		}

		@Override
		public String toString() {
			return "UmiReadInterval [CELL_BARCODE=" + CELL_BARCODE + ", GENE_NAME=" + GENE_NAME + ", MOLECULAR_BARCODE="
					+ MOLECULAR_BARCODE + ", NUM_READS=" + NUM_READS + ", NUM_SPLICED_READS=" + NUM_SPLICED_READS
					+ ", CONTIG=" + CONTIG + ", POSITION_MIN=" + POSITION_MIN + ", POSITION_MAX=" + POSITION_MAX
					+ ", POSITIVE_STRAND=" + POSITIVE_STRAND + "]";
		}		

		private String makeGapsString() {
			final String[] gapStrings = new String[GAPS.size()];
			int i = 0;
			for (Map.Entry<AlignmentGap, Integer> entry : GAPS.entrySet()) {
				AlignmentGap gap = entry.getKey();
				int count = entry.getValue();
				gapStrings[i++] = gap.start + ":" + gap.end + ":" + count;
			}
			return StringUtils.join(gapStrings, ",");
		}

		public String toTsv() {
			String [] body = {CELL_BARCODE, GENE_NAME, MOLECULAR_BARCODE, Integer.toString(NUM_READS), Integer.toString(NUM_SPLICED_READS),
					CONTIG, Integer.toString(POSITION_MIN), Integer.toString(POSITION_MAX), Boolean.toString(POSITIVE_STRAND),
					makeGapsString()};
			return StringUtils.join(body, "\t");
		}
	}
	
	
	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new GatherUMIReadIntervals().instanceMain(args));
	}
	
}
