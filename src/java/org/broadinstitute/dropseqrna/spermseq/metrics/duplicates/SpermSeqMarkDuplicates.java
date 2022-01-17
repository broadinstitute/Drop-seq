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
package org.broadinstitute.dropseqrna.spermseq.metrics.duplicates;

import htsjdk.samtools.*;
import htsjdk.samtools.DuplicateScoringStrategy.ScoringStrategy;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import org.apache.commons.math3.util.Precision;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.BarcodeListRetrieval;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.SpermSeq;
import org.broadinstitute.dropseqrna.utils.*;
import org.broadinstitute.dropseqrna.utils.readiterators.MissingTagFilteringIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamRecordSortingIteratorFactory;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.sam.DuplicationMetrics;

import java.io.File;
import java.util.*;

@CommandLineProgramProperties(
        summary = "For each cell/UMI, collect the set of reads." +
			"For reads that map to the same position and have the same cell and UMI, find the one with the best base quality and keep as the non duplicate record "+
			"Mark all other reads at that position/cell/UMI as duplicates.  This is similar to Picard's MarkDuplicates, but only works on single ended reads, and "
			+ "leverages UMIs to distinguish between reads mapped at the same position.",

        oneLineSummary = "Mark SpermSeq PCR Duplicates ",
        programGroup = SpermSeq.class
)
public class SpermSeqMarkDuplicates extends CommandLineProgram  {

	private static final Log log = Log.getInstance(SpermSeqMarkDuplicates.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze. This argument can accept wildcards, or a file with the suffix .bam_list that contains the locations of multiple BAM files", minElements = 1)
	public List<File> INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output SAM or BAM file with the PCRDuplicate flag set on reads.")
	public File OUTPUT;

	@Argument (doc="The file to output PCR duplication stats to", optional=false)
	public File OUTPUT_STATS;

	@Argument(doc="The cell barcode tag.")
	public String CELL_BARCODE_TAG="XC";

	@Argument(doc="The molecular barcode tag.")
	public String MOLECULAR_BARCODE_TAG="XM";

	@Argument(doc="The map quality of the read to be included when calculating whcih cells are the top cells to be reported.  "
			+ "This does not effect which reads are marked as duplicates, only which cells are selected for output")
	public Integer READ_MQ=10;

	@Argument (doc="Find the top set of cell barcodes by looking at the top <NUM_BARCODES> most common barcodes by HQ reads and only use this set for the reports.  If left unset, only the global summary is reported.", optional=true)
	public Integer NUM_BARCODES;

	@Argument (doc="Find the top set of cell barcodes by looking at all barcodes that have at least <NUM_READS> HQ reads. Only use this set for the reports.  If left unset, only the global summary is reported.", optional=true)
	public Integer NUM_READS;

	@Argument(doc="Override NUM_BARCODES and NUM_READS, and process reads that have the cell barcodes in this file instead.  The file has 1 column with no header.", optional=true)
	public File CELL_BC_FILE=null;

	@Argument (doc="Run with verbose output to help understand maximum memory usage and large pileups.")
	public boolean VERBOSE=false;

    @Argument(shortName = StandardOptionDefinitions.SORT_ORDER_SHORT_NAME)
    public SAMFileHeader.SortOrder SORT_ORDER = SAMFileHeader.SortOrder.coordinate;

	private final String POS_TAG="ZZ";
	private ScoringStrategy DUPLICATE_SCORING_STRATEGY = ScoringStrategy.SUM_OF_BASE_QUALITIES;
	private final String AGGREGATE_NAME="ALL";

	private PCRDuplicateMetrics globalMetrics=new PCRDuplicateMetrics();;

	private SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE);

	public enum DuplicateStrategy {
		READ_POSITION,
		CLUSTER;
	}
	
	
	@Override
	protected int doWork() {
		
		
		// for reporting, what is the aggregate output name
		globalMetrics.CELL_BARCODE=AGGREGATE_NAME;

		List<String> cellBarcodes = getCellBarcodes();

        // no need to maintain input sort because records must be re-sorted to do dupe marking
		SamHeaderAndIterator headerAndIterator = SamFileMergeUtil.mergeInputs(INPUT, false, samReaderFactory);
		SAMFileHeader h= headerAndIterator.header;
        h.setSortOrder(SORT_ORDER);
		SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(h, false, OUTPUT);

		Map<String, PCRDuplicateMetrics> metricsMap=null;

		metricsMap=processReadsByPosition(cellBarcodes, headerAndIterator, writer);

		log.info("Processing done, writing final BAM - this requires resorting the output BAM into genomic order");
		writer.close();

		MetricsFile<PCRDuplicateMetrics, String> file = new MetricsFile<>();
		this.globalMetrics.calculateStats();
		file.addMetric(this.globalMetrics);
		for (String key: cellBarcodes) {
			PCRDuplicateMetrics cell = metricsMap.get(key);
			if (cell!=null) {
				cell.calculateStats();
				file.addMetric(cell);
			}

		}
		file.write(this.OUTPUT_STATS);
		return 0;
	}

	/**
	 * Takes the merged SAM header + iterator, constructs a grouping iterator, then iterates through and markes each read as a duplicate or not duplicate read.
	 * Updates the metricsMap
	 * @param headerAndIterator
	 * @param writer
	 */
	Map<String, PCRDuplicateMetrics> processReadsByPosition (final List<String> cellBarcodes, final SamHeaderAndIterator headerAndIterator, final SAMFileWriter writer) {
		final Iterator<SAMRecord> filteringIterator = new MissingTagFilteringIterator(headerAndIterator.iterator, this.MOLECULAR_BARCODE_TAG, this.CELL_BARCODE_TAG);
	    // sort by Cell and molecular barcode
	    @SuppressWarnings("unchecked")
		final MultiComparator<SAMRecord> comparator = new MultiComparator<>(
	            new StringTagComparator(this.CELL_BARCODE_TAG), new StringTagComparator(this.MOLECULAR_BARCODE_TAG), new IntervalTagComparator(this.POS_TAG, headerAndIterator.header.getSequenceDictionary()));
	    // add the position tag.
	    final ReadDuplicateWrapper sortingIteratorWrapper = new ReadDuplicateWrapper(filteringIterator, POS_TAG);
	    final CloseableIterator<SAMRecord> sortingIterator =
	            SamRecordSortingIteratorFactory.create(headerAndIterator.header, sortingIteratorWrapper, comparator, new ProgressLogger(this.log));
	    final GroupingIterator<SAMRecord> groupingIterator = new GroupingIterator<>(sortingIterator, comparator);

		Map<String, PCRDuplicateMetrics> metricsMap = getPerCellMetricsMap(cellBarcodes);

		// logger to record how many records are being written out.
		ProgressLogger plWriter = new ProgressLogger(this.log, 1000000, "Wrote duplicate marked", "SAM records");

		int maxBatchSize=0;
		BatchStats stats=null;

		while (groupingIterator.hasNext()) {
			Collection<SAMRecord> batch = groupingIterator.next();
			if (batch.size()>maxBatchSize) {
				maxBatchSize=batch.size();
				stats=new BatchStats(this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, this.POS_TAG, batch.iterator().next(), batch.size());
				if (this.VERBOSE) log.info(stats);
			}
			Collection<SAMRecord> result = markDuplicates(batch, metricsMap);
			for (SAMRecord r: result) {
				r.setAttribute(this.POS_TAG, null);
				writer.addAlignment(r);
				plWriter.record(r);
			}
		}
		if (stats!=null) log.info("Maximum Batch Size " + stats);
		CloserUtil.close(groupingIterator);
		return metricsMap;
	}

	private class BatchStats {
		private int batchSize;
		private String cellBarcode;
		private String molBC;
		private String pos;

		public BatchStats (final String cellBarcodeTag, final String molBCTag, final String posTag, final SAMRecord r, final int batchSize) {
			this.batchSize=batchSize;
			this.cellBarcode = r.getStringAttribute(cellBarcodeTag);
			this.molBC = r.getStringAttribute(molBCTag);
			this.pos = r.getStringAttribute(posTag);
		}

		@Override
		public String toString () {
			return ("Cell Barcode [" + this.cellBarcode+ "] molBC ["+ this.molBC +"] position ["+ this.pos+"] batchSize ["+ this.batchSize+"]");
		}
	}

	/**
	 * Returns a map where the key is the cell barcode, and the value is the PCRDuplicateMetrics result for that barcode
	 * @param cellBarcodes
	 * @return
	 */
	Map<String, PCRDuplicateMetrics> getPerCellMetricsMap (final List<String> cellBarcodes) {
		Map<String, PCRDuplicateMetrics> map = new HashMap<>();
		for (String key: cellBarcodes) {
			PCRDuplicateMetrics m = new PCRDuplicateMetrics();
			m.CELL_BARCODE=key;
			map.put(key, m);
		}
		return (map);
	}

	/**
	 * Get a list of cell barcodes to report in the output.  If NUM_BARCODES is null, returns an empty list.
	 */
	List<String> getCellBarcodes () {

		if (this.CELL_BC_FILE!=null) {
			List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(this.CELL_BC_FILE);
			log.info("Found " + cellBarcodes.size()+ " cell barcodes in file");
			return (cellBarcodes);
		}

		if (this.NUM_BARCODES==null && this.NUM_READS==null) return new ArrayList<>();
		if (this.NUM_BARCODES!=null) {
			log.info("Gathering barcodes for the top [" + this.NUM_BARCODES +"] cells");

	        return new BarcodeListRetrieval().getListCellBarcodesByReadCount(
	                SamFileMergeUtil.mergeInputs(INPUT, false, samReaderFactory).iterator,
	                this.CELL_BARCODE_TAG, this.READ_MQ, null, this.NUM_BARCODES);
		}

		log.info("Gathering barcodes for cells with at least [" + this.NUM_READS +"] reads");
		return new BarcodeListRetrieval().getListCellBarcodesByReadCount(
	               	SamFileMergeUtil.mergeInputs(INPUT, false, samReaderFactory).iterator,
	                this.CELL_BARCODE_TAG, this.READ_MQ, this.NUM_READS, null);

	}





	/**
	 * Process a batch of reads
	 * In position strategy: that all have the same cell barcode, molecular barcode, and chromosome/untrimmed position.
	 * In cluster strategy: that all have the same cell barcode, molecular barcode and are within some distance from each other.
	 * Calculate the baseScore via the Picard method, and select the highest scoring read as not duplicated.
	 * @param batch A collection of reads to process
	 * @return The input read collection with the duplicate flag set appropriately on each read.
	 */
	Collection<SAMRecord> markDuplicates (final Collection<SAMRecord> batch, final Map<String, PCRDuplicateMetrics> metricsMap) {
		if (batch.size()==0) return (batch); // if there's 0 reads, you can't mark duplicates.

		// map read names to scores
		String topRead=null;
		int topScore=-1;
		// find max score, assign scores to reads.
		for (SAMRecord rec: batch) {

			int score = DuplicateScoringStrategy.computeDuplicateScore(rec, this.DUPLICATE_SCORING_STRATEGY);
			// log.info("Read Name [" + rec.getReadName()+ "] score [" + Integer.toString(score)+"]");
			// only look at mapped reads for the best read by base quality score.
			if (score>topScore && !rec.getReadUnmappedFlag()) {
				topRead=rec.getReadName();
				topScore=score;
			}
		}

		// assign the scores.
		for (SAMRecord rec: batch) {
			PCRDuplicateMetrics cell = metricsMap.get(rec.getAttribute(this.CELL_BARCODE_TAG));
			this.globalMetrics.NUM_READS++;
			if (cell!=null)
				cell.NUM_READS++;
			// only modify mapped reads.
			// only count mapped reads for duplicates.
			if (!rec.getReadUnmappedFlag()) {
				globalMetrics.NUM_MAPPED_READS++;
				if (cell!=null)
					cell.NUM_MAPPED_READS++;
				if (rec.getReadName().equals(topRead))
					rec.setDuplicateReadFlag(false);
				else {
					rec.setDuplicateReadFlag(true);
					globalMetrics.NUM_DUPLICATES++;
					if (cell!=null)
						cell.NUM_DUPLICATES++;
				}
			}
		}
		return batch;
	}





	public static class PCRDuplicateMetrics extends MetricBase {
		public String CELL_BARCODE;
		public int NUM_READS;
		public int NUM_MAPPED_READS;
		public int NUM_DUPLICATES;
		public double PCT_DUPLICATES;
        /** Note that this value is not particularly meaningful for summary (multi-library) metrics */
		public Long ESTIMATED_LIBRARY_SIZE;

		public void calculateStats() {
			this.PCT_DUPLICATES=(double) this.NUM_DUPLICATES / (double) this.NUM_MAPPED_READS * 100;
			this.PCT_DUPLICATES=Precision.round(this.PCT_DUPLICATES, 3);
			this.ESTIMATED_LIBRARY_SIZE = DuplicationMetrics.estimateLibrarySize(this.NUM_MAPPED_READS, this.NUM_MAPPED_READS - this.NUM_DUPLICATES);
		}
	}


	GroupingIterator<SAMRecord> getCellGeneIterator (final SamHeaderAndIterator headerAndIterator) {
		final Iterator<SAMRecord> filteringIterator = new MissingTagFilteringIterator(headerAndIterator.iterator, this.MOLECULAR_BARCODE_TAG, this.CELL_BARCODE_TAG);

	    // sort by Cell and molecular barcode
	    @SuppressWarnings("unchecked")
		final MultiComparator<SAMRecord> comparator = new MultiComparator<>(
	            new StringTagComparator(this.CELL_BARCODE_TAG), new StringTagComparator(this.MOLECULAR_BARCODE_TAG), new StringTagComparator(this.POS_TAG));

	    // add the position tag.
	    final ReadDuplicateWrapper sortingIteratorWrapper = new ReadDuplicateWrapper(filteringIterator, POS_TAG);

	    final CloseableIterator<SAMRecord> sortingIterator =
	            SamRecordSortingIteratorFactory.create(headerAndIterator.header, sortingIteratorWrapper, comparator, new ProgressLogger(this.log));

	    final GroupingIterator<SAMRecord> groupingIterator = new GroupingIterator<>(sortingIterator, comparator);
	    return groupingIterator;
	}
	
	private class ReadMappedFilteredIterator extends FilteredIterator<SAMRecord> {
		public ReadMappedFilteredIterator(final Iterator<SAMRecord> underlyingIterator) {
	        super(underlyingIterator);
		}
	    @Override
	    public boolean filterOut(final SAMRecord r) {
	        return r.getReadUnmappedFlag();
	    }
	}

	
	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new SpermSeqMarkDuplicates().instanceMain(args));
	}
	
	@Override
	protected String[] customCommandLineValidation() {

		final ArrayList<String> list = new ArrayList<>(1);
		this.INPUT = FileListParsingUtils.expandFileList(INPUT);
		
		for (File i: INPUT)
			IOUtil.assertFileIsReadable(i);
		
		IOUtil.assertFileIsWritable(OUTPUT);
		IOUtil.assertFileIsWritable(OUTPUT_STATS);
		if (this.CELL_BC_FILE!=null) 
			IOUtil.assertFileIsReadable(this.CELL_BC_FILE);
		
		return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);
		
	}



}
