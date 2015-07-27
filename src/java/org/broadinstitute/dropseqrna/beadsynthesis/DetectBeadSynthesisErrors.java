package org.broadinstitute.dropseqrna.beadsynthesis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
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
	
	@Option (doc="Find the top set of <NUM_BARCODES> most common barcodes by HQ reads and only use this set for analysis.")
	public Integer NUM_BARCODES;
		
	private Double EXTREME_BASE_RATIO=0.9;
	private Character PAD_CHARACTER='N';
	private Integer MAX_NUM_ERRORS=2;
	
	@Override
	protected int doWork() {
		
		
		//TagOrderIterator toi = checkInputsAndPrepIter();
		UMIIterator iter = checkInputsAndPrepIter();
		 
		// initialize output writers.
		// PrintStream outBarcodes = new PrintStream(IOUtil.openFileForWriting(OUT_BARCODES));
		PrintStream out = new PrintStream(IOUtil.openFileForWriting(OUTPUT_STATS));
		
		// for holding barcodes results.  The key is the barcode, the value is the first base to pad.
		// Used for cleanup of BAMs, if needed.
		Map<String, BeadSynthesisErrorData> errorBarcodesWithPositions = new HashMap<String, BeadSynthesisErrorData>();
		
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
			bsed.addUMI(umiCol);
			counter++;
			if (counter%1000000==0) log.info("Processed [" + counter + "] Cell/Gene UMIs.");
		}
		
		iter.close();
		
		// filter so these errors have a minimum number of UMIs.
		errorBarcodesWithPositions = filterByNumUMis(errorBarcodesWithPositions);
		
		writeFile (errorBarcodesWithPositions.values(), out);
		
		// clean up the BAM if desired.
		if (this.OUT!=null) {
			cleanBAM(errorBarcodesWithPositions, this.INPUT, this.OUT);
		}
		return 0;
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
			r=padCellBarcodeFix(r, errorBarcodesWithPositions, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG);
			if (r!=null) {
				writer.addAlignment(r);
			}
			
			pl.record(r);
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
		
		int errorPosition = bsed.getErrorBase(this.EXTREME_BASE_RATIO);
		int umiLength = bsed.getBaseLength();
		int numErrors= umiLength-errorPosition+1;
		if (numErrors > MAX_NUM_ERRORS) {
			return null;
		}
		String cellBCFixed = padCellBarcode(cellBC, errorPosition, umiLength);
		String umi = r.getStringAttribute(molecularBarcode);
		String umiFixed = fixUMI(cellBC, umi, errorPosition); 
		r.setAttribute(cellBarcodeTag, cellBCFixed);
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
	
	private void writeFile (Collection <BeadSynthesisErrorData> data, PrintStream out) {
		if (data.size()==0) return;
		
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
		header.add("error_base");
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
		//IOUtil.assertFileIsWritable(this.OUT_BARCODES);
		IOUtil.assertFileIsWritable(this.OUTPUT_STATS);
		
		if (OUT!=null) IOUtil.assertFileIsWritable(this.OUT);
		log.info("Gathering barcodes for the top [" + this.NUM_BARCODES +"] cells");
		List<String> barcodes = new BarcodeListRetrieval().getListCellBarcodesByReadCount (INPUT, this.CELL_BARCODE_TAG, this.READ_MQ, null, this.NUM_BARCODES);
		
		UMIIterator iter = new UMIIterator(this.INPUT, this.GENE_EXON_TAG, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, this.STRAND_TAG, this.READ_MQ, false, false, barcodes, this.MAX_RECORDS_IN_RAM);
		
		/*
		List<String> l = new ArrayList<String>();
		l.add(CELL_BARCODE_TAG);
		l.add(MOLECULAR_BARCODE_TAG);
		TagValueProcessor tvp = new TagValueProcessor(this.CELL_BARCODE_TAG, barcodes, true);
		ReadProcessorCollection rpc = new ReadProcessorCollection();
		rpc.addFilter(tvp);
		ComparatorAggregator ag = new ComparatorAggregator(new StringComparator(), true);
		TagOrderIterator toi = new TagOrderIterator(INPUT, l, l, ag, rpc, true);
		*/
		return (iter);
	}
	
	
	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new DetectBeadSynthesisErrors().instanceMain(args));
	}
}
