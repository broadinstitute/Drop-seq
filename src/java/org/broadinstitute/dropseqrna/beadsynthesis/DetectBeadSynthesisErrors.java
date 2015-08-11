package org.broadinstitute.dropseqrna.beadsynthesis;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.dropseqrna.barnyard.BarcodeListRetrieval;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.*;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

/**
 * 
 * @author nemesh
 *
 */

@CommandLineProgramProperties(
        usage = "For each cell, gather up all the UMIs.  An error in synthesis will result in the last base of the synthesis being fixed in >90% of the UMIs for that cell, across all genes." +
			"This fixed base is T.  For cell barcodes where this occurs, output the cell barcode in a file, as well as (optionally) pad the cell barcodes with N for the error bases.",
        usageShort = "Detect barcode synthesis errors where the final base of a UMI is fixed across all UMIs of a cell.",
        programGroup = DropSeq.class
)

public class DetectBeadSynthesisErrors extends CommandLineProgram {
	
	private static final Log log = Log.getInstance(DetectBeadSynthesisErrors.class);
	
	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM files to analyze.  They must all have the same sort order")
	public List<File> INPUT;
	
	@Option(doc="Output of detailed information on each cell barcode analyzed.  Each row is a single cell barcode.  "
			+ "The data has multiple columns: the cell barcode, the number of UMIs, then one column per UMI base position containing the count of the reads, with a | "
			+ "delimiter between bases.  Bases are ordered A,C,G,T for these columns.  An example output with a single base UMI would be:"
			+ "AAAAAA	20		5|4|6|5.")
	public File OUTPUT_STATS;
	
	@Option(doc="Output a summary of the error types and frequencies detected")
	public File SUMMARY;
	
	@Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output BAM, with the synthesis error barcodes removed", optional=true)
	public File OUTPUT;
	
	@Option(doc="The sequence of the primer.", optional=true)
	public String PRIMER_SEQUENCE=null;
	
	@Option(doc="When looking at fixed UMIs, see if the edit distance from the UMI to the primer is within this threshold.  0 indicates a perfect match between the primer and the UMI.")
	public Integer EDIT_DISTANCE=0;
	
	@Option(doc="The cell barcode tag.")
	public String CELL_BARCODE_TAG="XC";
	
	@Option(doc="The molecular barcode tag.")
	public String MOLECULAR_BARCODE_TAG="XM";
	
	@Option(doc="The Gene/Exon tag")
	public String GENE_EXON_TAG="GE";
	
	@Option(doc="The strand of the gene(s) the read overlaps.  When there are multiple genes, they will be comma-separated.")
	public String STRAND_TAG="GS";
	
	@Option(doc="The map quality of the read to be included when calculating the barcodes in <NUM_BARCODES>")
	public Integer READ_MQ=10;
	
	@Option (doc="The minimum number of UMIs required to report a cell barcode")
	public Integer MIN_UMIS_PER_CELL=25;
	
	@Option (doc="Find the top set of <NUM_BARCODES> most common barcodes by HQ reads and only use this set for analysis.",
            mutex = {"CELL_BC_FILE"})
	public Integer NUM_BARCODES;
	
	@Option(doc="Override NUM_BARCODES, and process reads that have the cell barcodes in this file instead.  The file has 1 column with no header.",
            mutex = {"NUM_BARCODES"})
	public File CELL_BC_FILE;
	
	@Option(doc="Repair Synthesis errors with at most this many missing bases detected.", optional=true)
	public Integer MAX_NUM_ERRORS=1;
	
	private Double EXTREME_BASE_RATIO=0.8;
	private Character PAD_CHARACTER='N';
	
	private DetectPrimerInUMI detectPrimerTool=null;

    private SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE);
	
	@Override
	protected int doWork() {
		// primer detection if requested.
		if (this.PRIMER_SEQUENCE!=null) {
			detectPrimerTool = new DetectPrimerInUMI(this.PRIMER_SEQUENCE);
		}

        // Group the stream of UMICollections into groups with the same cell barcode.
        GroupingIterator<UMICollection> groupingIterator = new GroupingIterator<>(prepareIterator(),
                new Comparator<UMICollection>() {
                    @Override
                    public int compare(UMICollection o1, UMICollection o2) {
                        return o1.getCellBarcode().compareTo(o2.getCellBarcode());
                    }
                });
		 

		// for holding barcodes results.  The key is the cell barcode, the value is the first base to pad.
		// Used for cleanup of BAMs, if needed.
		Map<String, BeadSynthesisErrorData> errorBarcodesWithPositions = new HashMap<String, BeadSynthesisErrorData>();

 		// for holding UMI Strings efficiently
		StringInterner  umiStringCache = new StringInterner();
		
		int counter=0;
        int numCellsFilteredLowUMIs = 0;

        for (final List<UMICollection> umiCollectionList : groupingIterator) {
            final String cellBarcode = umiCollectionList.get(0).getCellBarcode();
            BeadSynthesisErrorData bsed = new BeadSynthesisErrorData(cellBarcode);
            for (final UMICollection umis : umiCollectionList) {
                Collection<String> umiCol = umis.getMolecularBarcodes();
                umiCol=getUMIsFromCache(umiCol,umiStringCache);
                bsed.addUMI(umiCol);
                counter++;
                if (counter%1000000==0) log.info("Processed [" + counter + "] Cell/Gene UMIs.");
            }
            if (bsed.getUMICount() < this.MIN_UMIS_PER_CELL) {
                ++numCellsFilteredLowUMIs;
            } else {
                errorBarcodesWithPositions.put(cellBarcode, bsed);
            }
        }


		CloserUtil.close(groupingIterator);

        // initialize output writers.
        PrintStream out = new PrintStream(IOUtil.openFileForWriting(OUTPUT_STATS));
		writeFile (errorBarcodesWithPositions.values(), out);
		writeSummary(errorBarcodesWithPositions.values(), numCellsFilteredLowUMIs, SUMMARY);
		
		// clean up the BAM if desired.
		if (this.OUTPUT!=null) {
			cleanBAM(errorBarcodesWithPositions);
		}
		return 0;
	}
	
	/**
	 * Gets a reference to the UMI strings from the cache.  Has the side effect of populating the cache with additional 
	 * strings.  This reduces total memory footprint by returning references to repeated strings instead of 
	 * holding new objects for the same UMI over and over.
	 * @param umis A list of strings to get references to
	 * @param umiStringCache The cache of strings holding references.
	 */
	private Collection<String> getUMIsFromCache (Collection<String> umis, StringInterner umiStringCache) {
		List<String> result = new ArrayList<String>(umis.size());
		for (String umi: umis) {
            result.add(umiStringCache.intern(umi));
		}
		return (result);
	}
	
	/**
	 * For each problematic cell, replace cell barcodes positions with N.
	 * Take the replaced bases and prepend them to the UMI, and trim the last <X> bases off the end of the UMI.
	 */
	private void cleanBAM (Map<String, BeadSynthesisErrorData> errorBarcodesWithPositions) {
		log.info("Cleaning BAM");
        final SamHeaderAndIterator headerAndIterator = SamFileMergeUtil.mergeInputs(INPUT, true);
		SAMFileWriter writer= new SAMFileWriterFactory().setCreateIndex(CREATE_INDEX).makeSAMOrBAMWriter(headerAndIterator.header, true, OUTPUT);
		ProgressLogger pl = new ProgressLogger(log);
		for (SAMRecord r: new IterableAdapter<>(headerAndIterator.iterator)) {
			pl.record(r);
			r=padCellBarcodeFix(r, errorBarcodesWithPositions, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, this.EXTREME_BASE_RATIO);
			if (r!=null) {
				writer.addAlignment(r);
			}			
		}
		CloserUtil.close(headerAndIterator.iterator);
		writer.close();
	}
	
	
	/**
	 * @return null if the read should not be included in the output BAM.
	 */
	SAMRecord padCellBarcodeFix (SAMRecord r, Map<String, BeadSynthesisErrorData> errorBarcodesWithPositions, String cellBarcodeTag, String molecularBarcodeTag, double extremeBaseRatio) {
		String cellBC=r.getStringAttribute(cellBarcodeTag);
		
		BeadSynthesisErrorData bsed = errorBarcodesWithPositions.get(cellBC);
		if (bsed==null) return (r); // no correction data, no fix.
		
		// we're only going to fix cells where there's one or more synthesis errors
		BeadSynthesisErrorTypes bset = getEnhancedErrorType(bsed, extremeBaseRatio);
		if (bset==BeadSynthesisErrorTypes.NO_ERROR) return (r); // no error, return.
		// has an error, not a synthesis error...
		if (bset!=BeadSynthesisErrorTypes.SYNTH_MISSING_BASE) {
			return (null);
		}
		
		// has a synthesis error
		int polyTErrorPosition = bsed.getPolyTErrorPosition(this.EXTREME_BASE_RATIO);
		int umiLength = bsed.getBaseLength();
		int numErrors= umiLength-polyTErrorPosition+1;
		// if there are too many errors, or the errors aren't all polyT, return null.
		if (numErrors > MAX_NUM_ERRORS) {
			return null;			
		}
		
		// apply the fix and return the fixed read.
		String umi = r.getStringAttribute(molecularBarcodeTag);
		String cellBCFixed = padCellBarcode(cellBC, polyTErrorPosition, umiLength);
		String umiFixed = fixUMI(cellBC, umi, polyTErrorPosition); 
		r.setAttribute(cellBarcodeTag, cellBCFixed);
		r.setAttribute(molecularBarcodeTag, umiFixed);
		return r;
	}
	
	/**
	 * Take the original cell barcode and UMI, and move bases from the end of the cell barcode to the start of the UMI,
	 * then trim an equal number of bases off the end of the UMI so the length is the same.
	 * Example:
	 * Cell barcode: 		ACGCTCATACAG
	 * UMI: 				TCCTTATT
	 * errorPosition: 		2
	 * New Cell Barcode:	ACGCTCATACNN
	 * New UMI:				AGTCCTTA
	 * 
	 * @param cellBarcode The original cell barcode
	 * @param umi The original UMI 
	 * @param errorPosition The position in the UMI where the error occurred.
	 */
	String fixUMI (String cellBarcode, String umi, int errorPosition) {
		// 0 based, from end of cell barcode.
		int badBasesUMI=umi.length()-errorPosition;
		int lastBase = cellBarcode.length();
		int firstBaseToPad = lastBase-badBasesUMI-1;
		String cellBCBases=cellBarcode.substring(firstBaseToPad, cellBarcode.length());
		
		String umiRemaining=umi.substring(0, errorPosition-1);
		return cellBCBases+umiRemaining;
	}
	/**
	 * Picks a number of bases to pad.
	 * If errorPosition =-1, then don't pad any bases.
	 */
	String padCellBarcode (String cellBarcode, int errorPosition, int umiLength) {
		if (errorPosition==-1) return (cellBarcode);
		
		// 0 based, from end of cell barcode.
		int badBasesUMI=umiLength-errorPosition;
		int lastBase = cellBarcode.length();
		int firstBaseToPad = lastBase-badBasesUMI-1;
		
		char [] charAr = cellBarcode.toCharArray();
		for (int i=firstBaseToPad; i<lastBase; i++) {
			charAr[i]=this.PAD_CHARACTER;
		}
		return new String (charAr);
	}
	
	private void writeSummary (Collection <BeadSynthesisErrorData> data, int numCellsFilteredLowUMIs, File out) {
		// skip if no data.
		if (data.size()==0) {
			return;
		}
		
		// gather up error types.
		BeadSynthesisErrorsSummaryMetric m = new BeadSynthesisErrorsSummaryMetric();
		m.LOW_UMI_COUNT=numCellsFilteredLowUMIs;
		
		for (BeadSynthesisErrorData bsde: data) {
			BeadSynthesisErrorTypes t = getEnhancedErrorType(bsde, EXTREME_BASE_RATIO); 
			m.NUM_BEADS++;
			switch (t) {
				case SYNTH_MISSING_BASE: m.SYNTHESIS_MISSING_BASE++; m.incrementSynthesisMissingBase(bsde.getPolyTErrorPosition(this.EXTREME_BASE_RATIO)); break;
				case PRIMER: m.PRIMER_MATCH++; break;
				case SINGLE_UMI: m.SINGLE_UMI_ERROR++; break;
				case FIXED_FIRST_BASE: m.FIXED_FIRST_BASE++; break;
				case OTHER_ERROR: m.OTHER_ERROR_COUNT++; break;
				default: m.NO_ERROR++; break; 
			}
						
		}
		
		MetricsFile<BeadSynthesisErrorsSummaryMetric, Integer> outFile = new MetricsFile<BeadSynthesisErrorsSummaryMetric, Integer>();
		outFile.addMetric(m);
		outFile.addHistogram(m.getHistogram());
		outFile.write(out);
		
	}
	
	private void writeFile (Collection <BeadSynthesisErrorData> data, PrintStream out) {
		if (data.size()==0) {
			out.close();
			return;
		}
		
		// if there are records, write out the file.
		
		BeadSynthesisErrorData first = data.iterator().next();
		int umiLength = first.getBaseLength();
		writeBadBarcodeStatisticsFileHeader(umiLength, out);
		for (BeadSynthesisErrorData bsde: data) {
			writeBadBarcodeStatisticsFileEntry(bsde, out);
		}
		out.close();
	}
	
	/**
	 * Write the header.
	 */
	private void writeBadBarcodeStatisticsFileHeader (int umiLength, PrintStream out) {
		List<String> header = new ArrayList<String>();
		header.add("CELL_BARCODE");
		header.add("NUM_UMI");
		header.add("FIRST_BIASED_BASE");
		header.add(BeadSynthesisErrorTypes.SYNTH_MISSING_BASE.toString());
		header.add("ERROR_TYPE");
		for (int i=0; i<umiLength; i++) {
			header.add("BASE_"+ Integer.toString(i+1));
		}
		String h = StringUtils.join(header, "\t");
		out.println(h);
	}
	
	private void writeBadBarcodeStatisticsFileEntry (BeadSynthesisErrorData data, PrintStream out) {
		List<String> line = new ArrayList<String>();
		line.add(data.getCellBarcode());
		line.add(Integer.toString(data.getUMICount()));
		int base = data.getErrorBase(EXTREME_BASE_RATIO);
		line.add(Integer.toString(base));
		int polyTErrorBase = data.getPolyTErrorPosition(EXTREME_BASE_RATIO);
		line.add(Integer.toString(polyTErrorBase));
		line.add(getEnhancedErrorType(data, EXTREME_BASE_RATIO).toString());
		
		
		BaseDistributionMetricCollection bases = data.getBaseCounts();
		List<Integer> pos =bases.getPositions();
		for (Integer i: pos) {
			BaseDistributionMetric bdm = bases.getDistributionAtPosition(i);
			String formattedResult = format(bdm);
			line.add(formattedResult);
		}
		String outLine = StringUtils.join(line, "\t");
		out.println(outLine);
	}
	
	private BeadSynthesisErrorTypes getEnhancedErrorType (BeadSynthesisErrorData data, double extremeBaseRatio) {
		BeadSynthesisErrorTypes errorType = data.getErrorType(extremeBaseRatio);
		//base case, error is not a single UMI.
		if (errorType!=BeadSynthesisErrorTypes.SINGLE_UMI) {
			return errorType;
		}
		// if there's a primer, run detection.
		if (this.detectPrimerTool!=null) {
			// a single UMI-style error, does the most common UMI match the primer?
			String singleUMI = data.getUMICounts().getKeysOrderedByCount(true).get(0);
			boolean primerDetected = detectPrimerTool.isStringInPrimer(singleUMI, this.EDIT_DISTANCE);
			if (primerDetected) {
				return BeadSynthesisErrorTypes.PRIMER;
			} else {
				return errorType;
			}
		}
		
		 
		return errorType;
	}
	
	private String format (BaseDistributionMetric bdm) {
		List<String> d = new ArrayList<String>();
		
		for (Bases b: Bases.values()) {
			char bb = b.getBase();
			int count = bdm.getCount(bb);
			d.add(Integer.toString(count));
		}
		return StringUtils.join(d, "|");
	}

    @Override
    protected String[] customCommandLineValidation() {
        for (final File input : INPUT) {
            IOUtil.assertFileIsReadable(input);
        }
        IOUtil.assertFileIsWritable(this.OUTPUT_STATS);
        IOUtil.assertFileIsWritable(this.SUMMARY);

        if (OUTPUT!=null) IOUtil.assertFileIsWritable(this.OUTPUT);
        return super.customCommandLineValidation();
    }

    private UMIIterator prepareIterator() {
		List<String> barcodes=getCellBarcodes();
		return new UMIIterator(SamFileMergeUtil.mergeInputs(INPUT, false, samReaderFactory),
                this.GENE_EXON_TAG, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, this.STRAND_TAG, this.READ_MQ,
                false, false, barcodes, true);
	}
	
	private List<String> getCellBarcodes () {
		
		if (this.CELL_BC_FILE!=null) {
			IOUtil.assertFileIsReadable(this.CELL_BC_FILE);
			List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(this.CELL_BC_FILE);
			log.info("Found " + cellBarcodes.size()+ " cell barcodes in file");
			return (cellBarcodes);
		} 
		log.info("Gathering barcodes for the top [" + this.NUM_BARCODES +"] cells");
        return new BarcodeListRetrieval().getListCellBarcodesByReadCount(
                SamFileMergeUtil.mergeInputs(INPUT, false, samReaderFactory).iterator,
                this.CELL_BARCODE_TAG, this.READ_MQ, null, this.NUM_BARCODES);
	}
	
	public class BeadSynthesisErrorsSummaryMetric extends MetricBase {
		public int NUM_BEADS;
		public int NO_ERROR;
		public int LOW_UMI_COUNT;
		public int SYNTHESIS_MISSING_BASE;
		public int SINGLE_UMI_ERROR;
		public int PRIMER_MATCH;
		public int FIXED_FIRST_BASE;
		public int OTHER_ERROR_COUNT;
		
		
		/** The distribution of  SYNTHESIS_MISSING_BASE error positions */
		private Histogram <Integer> histogram = null;
				
		public BeadSynthesisErrorsSummaryMetric () {
			this.NUM_BEADS=0;
			this.NO_ERROR=0;
			this.SYNTHESIS_MISSING_BASE=0;
			this.SINGLE_UMI_ERROR=0;
			this.PRIMER_MATCH=0;
			this.OTHER_ERROR_COUNT=0;			
			this.LOW_UMI_COUNT=0;
			histogram = new Histogram<Integer>("SYNTHESIS_ERROR_BASE", "num cells");
		}
		
		public void incrementSynthesisMissingBase (int position) {
			histogram.increment(position);
		}
		
		public Histogram<Integer> getHistogram() {
			return histogram;
		}
			
	}
	
	
	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new DetectBeadSynthesisErrors().instanceMain(args));
	}
}
