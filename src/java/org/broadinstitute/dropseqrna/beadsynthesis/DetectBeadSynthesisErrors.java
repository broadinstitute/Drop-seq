package org.broadinstitute.dropseqrna.beadsynthesis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.dropseqrna.barnyard.BarcodeListRetrieval;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.BaseDistributionMetric;
import org.broadinstitute.dropseqrna.utils.BaseDistributionMetricCollection;
import org.broadinstitute.dropseqrna.utils.Bases;
import org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

/**
 * 
 * @author nemesh
 *
 */

@CommandLineProgramProperties(
        usage = "For each cell, gather up all the UMIs.  An error in synthesis will result in the last base of the synthesis being fixed in >90% of the UMIs for that cell, across all genes." +
			"This fixed base is T.  For cell barcodes where this occurs, output the cell barcode in a file, as well as (optinally) bad the cell barcodes with N for the error bases.", 
        usageShort = "Detect barconde synthesis errors where the final base of a UMI is fixed across all UMIs of a cell.",
        programGroup = DropSeq.class
)

public class DetectBeadSynthesisErrors extends CommandLineProgram {
	
	private static final Log log = Log.getInstance(DetectBeadSynthesisErrors.class);
	
	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.")
	public File INPUT;
	
	//@Option(doc="A list of all barcodes that are flagged as synthesis errors.  File has 1 column and no header.")
	//public File OUT_BARCODES;
	
	@Option(doc="Output of detailed information on each cell barcode analyzed.  Each row is a single cell barcode.  "
			+ "The data has multiple columns: the cell barcode, the number of umis, then one column per UMI base position containing the count of the reads, with a | "
			+ "delimiter between bases.  Bases are ordered A,C,G,T for these columns.  An example output with a single base UMI would be:"
			+ "AAAAAA	20		5|4|6|5.")
	public File OUTPUT_STATS;
	
	@Option(doc="Output a summary of the error types and frequencies detected")
	public File SUMMARY;
	
	@Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output BAM, with the synthesis error barcodes removed", optional=true)
	public File OUT;
	
	@Option(doc="The cell barcode tag.")
	public String CELL_BARCODE_TAG="XC";
	
	@Option(doc="The molecular barcode tag.")
	public String MOLECULAR_BARCODE_TAG="XM";
	
	@Option(doc="The Gene/Exon tag")
	public String GENE_EXON_TAG="GE";
	
	@Option(doc="The strand of the gene(s) the read overlaps.  When there are multiple genes, they will be comma seperated.")
	public String STRAND_TAG="GS";
	
	@Option(doc="The map quality of the read to be included when calculating the barcodes in <NUM_BARCODES>")
	public Integer READ_MQ=10;
	
	@Option (doc="The minimum number of UMIs required to report a cell barcode")
	public Integer MIN_UMIS_PER_CELL=25;
	
	@Option (doc="Find the top set of <NUM_BARCODES> most common barcodes by HQ reads and only use this set for analysis.", optional=true)
	public Integer NUM_BARCODES;
	
	@Option(doc="Override NUM_BARCODES, and process reads that have the cell barcodes in this file instead.  The file has 1 column with no header.", optional=true)
	public File CELL_BC_FILE=null;
		
	private Double EXTREME_BASE_RATIO=0.8;
	private Character PAD_CHARACTER='N';
	private Integer MAX_NUM_ERRORS=2;
	
	@Override
	protected int doWork() {
		// one or the other...
		if (this.CELL_BC_FILE==null & this.NUM_BARCODES==null) {
			log.error("Must set either CELL_BC_FILE or NUM_BARCODES");
			return 1;
		}		
		
		//TagOrderIterator toi = checkInputsAndPrepIter();
		UMIIterator iter = checkInputsAndPrepIter();
		 
		// initialize output writers.
		// PrintStream outBarcodes = new PrintStream(IOUtil.openFileForWriting(OUT_BARCODES));
		PrintStream out = new PrintStream(IOUtil.openFileForWriting(OUTPUT_STATS));
		
		// for holding barcodes results.  The key is the barcode, the value is the first base to pad.
		// Used for cleanup of BAMs, if needed.
		Map<String, BeadSynthesisErrorData> errorBarcodesWithPositions = new HashMap<String, BeadSynthesisErrorData>();
		// for holding UMI Strings
		Map<String, String> umiStringCache = new HashMap<String, String> ();
		
		int counter=0;
		
		while (iter.hasNext()) {
			UMICollection umis = iter.next();
			String cellBC = umis.getCellBarcode();
			
			Collection<String> umiCol = umis.getMolecularBarcodes();
			
			BeadSynthesisErrorData bsed = errorBarcodesWithPositions.get(cellBC);
			if (bsed==null) {
				bsed = new BeadSynthesisErrorData(cellBC);
				errorBarcodesWithPositions.put(cellBC, bsed);
			}
			// reduce memory BeadSynthesisErrorData objects by references to strings.
			umiCol=getUMIsFromCache (umiCol,umiStringCache);
			bsed.addUMI(umiCol);
			counter++;
			if (counter%1000000==0) log.info("Processed [" + counter + "] Cell/Gene UMIs.");
		}
		
		iter.close();
		
		// filter so these errors have a minimum number of UMIs.
		errorBarcodesWithPositions = filterByNumUMis(errorBarcodesWithPositions);
		
		writeFile (errorBarcodesWithPositions.values(), out);
		writeSummary(errorBarcodesWithPositions.values(), SUMMARY);
		
		// clean up the BAM if desired.
		if (this.OUT!=null) {
			cleanBAM(errorBarcodesWithPositions, this.INPUT, this.OUT);
		}
		return 0;
	}
	
	/**
	 * Gets a reference to the UMI strings from the cache.  Has the side effect of populating the cache with additional 
	 * strings.  This reduces total memory footprint by returning references to repeated strings instead of 
	 * holding new objects for the same UMI over and over.
	 * @param umis A list of strings to get references to
	 * @param umiStringCache The cache of strings holding references.
	 * @return
	 */
	private Collection<String> getUMIsFromCache (Collection<String> umis, Map<String, String> umiStringCache) {
		List<String> result = new ArrayList<String>(umis.size());
		for (String umi: umis) {
			String r = umiStringCache.get(umi);
			if (r==null) {
				umiStringCache.put(umi, umi);
				r=umi;
			}
			result.add(r);
		}
		return (result);
	}
	
	private Map<String, BeadSynthesisErrorData> filterByNumUMis (Map<String, BeadSynthesisErrorData> errorBarcodesWithPositions) {
		Map<String, BeadSynthesisErrorData> result = new HashMap<String, BeadSynthesisErrorData>();
		
		for (String cellBC: errorBarcodesWithPositions.keySet()) {
			BeadSynthesisErrorData bsed = errorBarcodesWithPositions.get(cellBC);
			if (bsed.getUMICount()>=this.MIN_UMIS_PER_CELL) {
				result.put(cellBC, bsed);
			}
		}
		
		return (result);
	}
	
	
	
	/**
	 * For each problematic cell, replace cell barcodes positions with N.
	 * Take the replaced bases and prepend them to the UMI, and trim the last <X> bases off the end of the UMI.
	 * @param errorBarcodesWithPositions
	 * @param inBAM
	 * @param outBAM
	 */
	private void cleanBAM (Map<String, BeadSynthesisErrorData> errorBarcodesWithPositions, File inBAM, File outBAM) {
		log.info("Cleaning BAM");
		SamReader reader = SamReaderFactory.makeDefault().open(inBAM);
		SAMFileHeader h= reader.getFileHeader();
		SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(h, true, outBAM);
		ProgressLogger pl = new ProgressLogger(this.log);
		for (SAMRecord r: reader) {
			pl.record(r);
			r=padCellBarcodeFix(r, errorBarcodesWithPositions, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG);
			if (r!=null) {
				writer.addAlignment(r);
			}			
		}
		CloserUtil.close(reader);
		CloserUtil.close(writer);
	}
	
	
	/**
	 * Returns null if the read should not be included in the output BAM.
	 * @param r
	 * @param errorBarcodesWithPositions
	 * @param cellBarcodeTag
	 * @param molecularBarcode
	 * @return
	 */
	SAMRecord padCellBarcodeFix (SAMRecord r, Map<String, BeadSynthesisErrorData> errorBarcodesWithPositions, String cellBarcodeTag, String molecularBarcode) {
		String cellBC=r.getStringAttribute(cellBarcodeTag);
		BeadSynthesisErrorData bsed = errorBarcodesWithPositions.get(cellBC);
		if (bsed==null) return (r); // no correction data, no fix.
		
		// we're only going to fix cells where the error position and the polyA error agree.
		// Essentially, we expect to fix 1 or 2 polyT's at the end of UMIs,
		// but other error modes should just be reported.
		// they will have separate cleanup methods if we can figure out how to correct them.
		int errorPosition = bsed.getErrorBase(this.EXTREME_BASE_RATIO);
		int polyAErrorPosition = bsed.getPolyTErrorPosition(this.EXTREME_BASE_RATIO);
		int umiLength = bsed.getBaseLength();
		int numErrors= umiLength-polyAErrorPosition+1;
		// if there are too many errors, or the errors aren't all polyT, return null.
		if (errorPosition!=polyAErrorPosition || numErrors > MAX_NUM_ERRORS) {
			return null;			
		}
		
		// apply the fix and return the fixed read.
		String umi = r.getStringAttribute(molecularBarcode);
		String cellBCFixed = padCellBarcode(cellBC, errorPosition, umiLength);
		String umiFixed = fixUMI(cellBC, umi, errorPosition); 
		r.setAttribute(cellBarcodeTag, cellBCFixed);
		r.setAttribute(molecularBarcode, umiFixed);
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
	 * @return
	 */
	String fixUMI (String cellBarcode, String umi, int errorPosition) {
		// 0 based, from end of cell barcode.
		int badBasesUMI=umi.length()-errorPosition;
		int lastBase = cellBarcode.length();
		int firstBaseToPad = lastBase-badBasesUMI-1;
		String cellBCBases=cellBarcode.substring(firstBaseToPad, cellBarcode.length());
		
		
		String umiRemaining=umi.substring(0, errorPosition-1);
		String newUMI=cellBCBases+umiRemaining;
		return (newUMI);
	}
	/**
	 * Picks a number of bases to pad.
	 * If errorPosition =-1, then don't pad any bases.
	 * @param cellBarcode
	 * @param errorPosition
	 * @param umiLength
	 * @return
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
		String fixedCellBarcode = new String (charAr);
		return (fixedCellBarcode);
	}
	
	private void writeSummary (Collection <BeadSynthesisErrorData> data, File out) {
		// skip if no data.
		if (data.size()==0) {
			return;
		}
		
		// gather up error types.
		BeadSynthesisErrorsSummaryMetric m = new BeadSynthesisErrorsSummaryMetric();
		
		for (BeadSynthesisErrorData bsde: data) {
			String eet = getEnhancedErrorType(bsde, EXTREME_BASE_RATIO);			
			BeadSynthesisErrorTypes et = bsde.getErrorType(this.EXTREME_BASE_RATIO);
			m.NUM_BEADS++;
			if (eet.equals(BeadSynthesisErrorTypes.POLY_T_ERROR.toString())) {
				m.POLY_T_ERROR_COUNT++;
			}
			else if (eet.equals(BeadSynthesisErrorTypes.SINGLE_UMI.toString())) {
				m.SINGLE_UMI_ERROR++; 
			}
			else if (eet.equals(BeadSynthesisErrorTypes.OTHER_ERROR.toString())) {
				m.OTHER_ERROR_COUNT++; 
			} 
			else if (eet.equals("PRIMER_MATCH")) {
				m.PRIMER_MATCH++;
			}
			else {
				m.NO_ERROR++;
			}			
		}
		
		MetricsFile<BeadSynthesisErrorsSummaryMetric, Integer> outFile = new MetricsFile<BeadSynthesisErrorsSummaryMetric, Integer>();
		outFile.addMetric(m);
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
	 * @param umiLength
	 * @param out
	 */
	private void writeBadBarcodeStatisticsFileHeader (int umiLength, PrintStream out) {
		List<String> header = new ArrayList<String>();
		header.add("cell_barcode");
		header.add("num_umi");
		header.add("any_error_base");
		header.add("polyT_error_base");
		header.add("error_type");
		for (int i=0; i<umiLength; i++) {
			header.add("base_"+ Integer.toString(i+1));
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
		line.add(getEnhancedErrorType(data, EXTREME_BASE_RATIO));
		
		
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
	
	private String getEnhancedErrorType (BeadSynthesisErrorData data, double extremeBaseRatio) {
		String errorType = data.getErrorType(extremeBaseRatio).toString();
		
		return errorType;
	}
	
	private String format (BaseDistributionMetric bdm) {
		List<String> d = new ArrayList<String>();
		
		for (Bases b: Bases.values()) {
			char bb = b.getBase();
			int count = bdm.getCount(bb);
			d.add(Integer.toString(count));
		}
		String result = StringUtils.join(d, "|");
		return result;
		
	}
	
	private UMIIterator checkInputsAndPrepIter () {
		IOUtil.assertFileIsReadable(this.INPUT);
		IOUtil.assertFileIsWritable(this.OUTPUT_STATS);
		IOUtil.assertFileIsWritable(this.SUMMARY);
		
		if (OUT!=null) IOUtil.assertFileIsWritable(this.OUT);
		List<String> barcodes=getCellBarcodes();		
		UMIIterator iter = new UMIIterator(this.INPUT, this.GENE_EXON_TAG, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, this.STRAND_TAG, this.READ_MQ, false, false, barcodes, this.MAX_RECORDS_IN_RAM);		
		return (iter);
	}
	
	private List<String> getCellBarcodes () {
		
		if (this.CELL_BC_FILE!=null) {
			IOUtil.assertFileIsReadable(this.CELL_BC_FILE);
			List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(this.CELL_BC_FILE);
			log.info("Found " + cellBarcodes.size()+ " cell barcodes in file");
			return (cellBarcodes);
		} 
		log.info("Gathering barcodes for the top [" + this.NUM_BARCODES +"] cells");
		List<String> barcodes = new BarcodeListRetrieval().getListCellBarcodesByReadCount (INPUT, this.CELL_BARCODE_TAG, this.READ_MQ, null, this.NUM_BARCODES);
		return(barcodes);
	}
	
	public class BeadSynthesisErrorsSummaryMetric extends MetricBase {
		public int NUM_BEADS;
		public int NO_ERROR;
		public int POLY_T_ERROR_COUNT;
		public int SINGLE_UMI_ERROR;
		public int PRIMER_MATCH;
		public int OTHER_ERROR_COUNT;
		
		public BeadSynthesisErrorsSummaryMetric () {
			this.NUM_BEADS=0;
			this.NO_ERROR=0;
			this.POLY_T_ERROR_COUNT=0;
			this.SINGLE_UMI_ERROR=0;
			this.OTHER_ERROR_COUNT=0;
		}
			
	}
	
	
	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new DetectBeadSynthesisErrors().instanceMain(args));
	}
}
