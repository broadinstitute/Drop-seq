package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.DigitalExpression;
import org.broadinstitute.dropseqrna.barnyard.GeneFunctionCommandLineBase;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        summary = "Filter BAM by the number of reads supporting a UMI."
        		+ "This preprocessing step can remove UMIs with less than or greater than some number of supporting reads.",                
        oneLineSummary = "Filter BAM by the number of reads supporting a UMI",
        programGroup = DropSeq.class
)

public class FilterReadsByUMISupport extends GeneFunctionCommandLineBase {

    private static final Log log = Log.getInstance(DigitalExpression.class);
    
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze. This argument can accept wildcards, or a file with the suffix .bam_list that contains the locations of multiple BAM files", minElements = 1)
	public List<File> INPUT;
	
	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output BAM File")
	public File OUTPUT;
	    
    @Argument(doc="A summary of the number of UMIs retained and removed", optional=true)
    public File METRICS;
    
    @Argument(shortName = StandardOptionDefinitions.SORT_ORDER_SHORT_NAME, doc="The sort order of the output BAM.  Many drop-seq tools that look at UMIs will sort the data, so it is not neccesary" +
    "to presort the BAM.  Defaults to unsorted, which saves time.")
    public SAMFileHeader.SortOrder SORT_ORDER = SAMFileHeader.SortOrder.unsorted;
    
    @Argument(doc="The minimum number of reads a UMI is supported by to be retained.  Set to -1 to ignore this filter")
    public Integer MIN_READ_SUPPORT=-1;
    
    @Argument(doc="The minimum number of reads a UMI is supported by to be retained.  Set to -1 to ignore this filter")
    public Integer MAX_READ_SUPPORT=-1;
    
    @Argument(doc="The cell barcode tag.  If there are no reads with this tag, the program will assume that all reads belong to the same cell and process in single sample mode.")
	public String CELL_BARCODE_TAG="XC";
	
	@Argument(doc="The molecular barcode tag.")
	public String MOLECULAR_BARCODE_TAG="XM";
	
	@Argument(doc="The edit distance that molecular barcodes should be combined at within a gene.")
	public Integer EDIT_DISTANCE=1;
	
	@Argument(doc="The map quality of the read to be included.")
	public Integer READ_MQ=10;
		
	@Argument(doc="If provided, process reads that have the cell barcodes in this file.  If not provided, all cells are processed. The file has 1 column with no header.", optional=true)
	public File CELL_BC_FILE=null;
    
    @Override
	protected int doWork() {
    	
    	List<String> cellBarcodes = getCellBarcodes();
    	    	
    	SamHeaderAndIterator headerAndIter = SamFileMergeUtil.mergeInputs(this.INPUT, false, SamReaderFactory.makeDefault());
    	SAMFileWriter out = getSamWriter (headerAndIter.header, this.SORT_ORDER);
    	
    	FilteredUmiMetrics metrics = new FilteredUmiMetrics();
    	    	
    	UMIIterator umiIterator = new UMIIterator.UMIIteratorBuilder(headerAndIter,GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG,
    			this.STRAND_STRATEGY, this.LOCUS_FUNCTION_LIST, this.FUNCTIONAL_STRATEGY, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG,
    			this.READ_MQ).setCellBarcodes(cellBarcodes).retainReads(true).build();
        
    	while (umiIterator.hasNext()) {
    		UMICollection c = umiIterator.next();
    		ObjectCounter<String> molecularBarcodeCounts = c.getMolecularBarcodeCountsCollapsed(this.EDIT_DISTANCE);
    		for (String molecularBarcode: molecularBarcodeCounts.getKeys()) {
    			metrics.NUM_UMIS++;
    			int count = molecularBarcodeCounts.getCountForKey(molecularBarcode);
    			if (retainCount (count, this.MIN_READ_SUPPORT, this.MAX_READ_SUPPORT)) {
    				// UMI retained, dump reads to bam.
    				List<SAMRecord> recs = c.getReads(molecularBarcode);
    				recs.forEach(x -> out.addAlignment(x));    				
    			} else {  
    				// UMI rejected
    				metrics.NUM_UMIS_REMOVED++;
    				metrics.histogram.increment(count);
    			}
    		}
    	}
        
    	log.info("Total cell/UMIs observed [" +metrics.NUM_UMIS+"] UMIs filtered [" + metrics.NUM_UMIS_REMOVED+"] % ["+ String.format("%.2f%%",metrics.getFractionMarked()*100)+"]");
    	CloserUtil.close(umiIterator);
    	CloserUtil.close(out);
    	
    	if (METRICS != null) {
            final MetricsFile<FilteredUmiMetrics, Integer> metricsFile = new MetricsFile<>();
            metricsFile.addMetric(metrics);
            metricsFile.addHistogram(metrics.getHistogram());
            metricsFile.write(METRICS);
        }
		return 0;
	}
    
    /**
     * Get a SAM writer.  If the output is unsorted, we can assume the data is presorted as it's added to the output writer.
     * @param header The input header to copy into the output file
     * @return A SAM writer in the requested sort order.
     */
    private SAMFileWriter getSamWriter (SAMFileHeader header, SAMFileHeader.SortOrder order) {
    	// is the data presorted?  If the output is unsorted, then yeah, kinda!
    	boolean presorted=order==SAMFileHeader.SortOrder.unsorted;    	

    	// set the sort order to the desired output order.
    	SAMFileHeader outHeader= header.clone();    	
    	outHeader.setSortOrder(this.SORT_ORDER);    	
    	SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(outHeader, presorted, OUTPUT);
    	return out;    	    	    	
    }
    
    /**
     * Tests if a count is in-between a lower and upper bound
     * @param count The count to test
     * @param lowerBound The lower bound. Count must be >= this value
     * @param upperBound The upper bound. Count must be <= this value
     * @return Returns true if the UMI number of supporting reads is between the lower and upper bounds
     */
    private boolean retainCount (int count, int lowerBound, int upperBound) {
    	return (lowerBound==-1 || count >= lowerBound) & (upperBound==-1 || count <= upperBound);    	
    }
    
    /**
     * Optionally parse the cell barcode list if not null.
     * @return A list of cell barcodes if a file is provided, or null if no list was provided. 
     */
    private List<String> getCellBarcodes () {
    	if (this.CELL_BC_FILE==null)
    		return null;
        
    	List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(CELL_BC_FILE);    	
    	log.info("Found " + cellBarcodes.size()+ " cell barcodes in file");
    	return cellBarcodes;    	        	
    }
    
    public static class FilteredUmiMetrics extends MetricBase {
		public long NUM_UMIS=0;
		public long NUM_UMIS_REMOVED=0;
		public double FRAC_UMIS_REMOVED;
		
		private Histogram <Integer> histogram = null;
		
		public FilteredUmiMetrics () {
			histogram = new Histogram<>("UMI_READ_COUNT", "NUM_UMIS_FILTERED");
		}
		
		public double getFractionMarked () {
			FRAC_UMIS_REMOVED= ((double) NUM_UMIS_REMOVED / (double) NUM_UMIS);
			return FRAC_UMIS_REMOVED;
		}

		/**
		 * Adds the non-calculated metrics (which is all of them)
		 */
		public void merge(final FilteredUmiMetrics metric) {
			this.NUM_UMIS += metric.NUM_UMIS;
			this.NUM_UMIS_REMOVED += metric.NUM_UMIS_REMOVED;	
			this.histogram.addHistogram(metric.histogram);
		}
		
		public Histogram<Integer> getHistogram() {
			return histogram;
		}
	}
    
    @Override
	protected String[] customCommandLineValidation() {
		final ArrayList<String> list = new ArrayList<>(1);		
		
		if (this.METRICS!=null) IOUtil.assertFileIsWritable(this.METRICS);
		if (this.CELL_BC_FILE!=null) IOUtil.assertFileIsReadable(this.CELL_BC_FILE);
		
		this.INPUT = FileListParsingUtils.expandFileList(INPUT);
		IOUtil.assertFileIsWritable(this.OUTPUT);
		
		if (this.MIN_READ_SUPPORT!=-1) {
			if (this.MIN_READ_SUPPORT<0) {
				list.add("If MIN_READ_SUPPORT is set, the value must be >=0.");
			}
		}
		
		if (this.MAX_READ_SUPPORT!=-1) {
			if (this.MAX_READ_SUPPORT<1) {
				list.add("If MAX_READ_SUPPORT is set, the value must be >=1.");
			}
		}					
		return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);
	}
    
    /** Stock main method. */
    public static void main(final String[] args) {
        System.exit(new FilterReadsByUMISupport().instanceMain(args));
    }

	
}


