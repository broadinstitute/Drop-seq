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
package org.broadinstitute.dropseqrna.beadsynthesis;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.broadinstitute.dropseqrna.barnyard.BarcodeListRetrieval;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.StringInterner;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator;

import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

public abstract class AbstractDetectBeadSynthesisErrors extends CommandLineProgram {

	private static final Log log = Log.getInstance(AbstractDetectBeadSynthesisErrors.class);

	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM files to analyze.  They must all have the same sort order")
	public List<File> INPUT;

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

	SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE);

	private Character PAD_CHARACTER='N';

	/**
	 * Find all the cell barcodes that are biased.
	 * @param iter
	 * @return
	 */
	public BiasedBarcodeCollection findBiasedBarcodes (final UMIIterator iter) {

		// Group the stream of UMICollections into groups with the same cell barcode.
        GroupingIterator<UMICollection> groupingIterator = new GroupingIterator<>(iter,
                new Comparator<UMICollection>() {
                    @Override
                    public int compare(final UMICollection o1, final UMICollection o2) {
                        return o1.getCellBarcode().compareTo(o2.getCellBarcode());
                    }
                });


		// for holding barcodes results.  The key is the cell barcode, the value is the first base to pad.
		// Used for cleanup of BAMs, if needed.
		Map<String, BeadSynthesisErrorData> errorBarcodesWithPositions = new HashMap<>();

 		// for holding UMI Strings efficiently
		StringInterner  umiStringCache = new StringInterner();

		int counter=0;
        int numCellsFilteredLowUMIs = 0;

        for (final List<UMICollection> umiCollectionList : groupingIterator) {
            final String cellBarcode = umiCollectionList.get(0).getCellBarcode();
            BeadSynthesisErrorData bsed = new BeadSynthesisErrorData(cellBarcode);
            for (final UMICollection umis : umiCollectionList) {
            	int transcriptCounts= umis.getDigitalExpression(1, 1, false);
            	int readCounts = umis.getDigitalExpression(1, 1, true);
                Collection<String> umiCol = umis.getMolecularBarcodes();
                umiCol=getUMIsFromCache(umiCol,umiStringCache);
                bsed.addUMI(umiCol);
                bsed.incrementReads(readCounts);
                bsed.incrementTranscripts(transcriptCounts);
                counter++;
                if (counter%1000000==0) log.info("Processed [" + counter + "] Cell/Gene UMIs.");
            }
            if (bsed.getUMICount() < this.MIN_UMIS_PER_CELL)
				++numCellsFilteredLowUMIs;
			else
				errorBarcodesWithPositions.put(cellBarcode, bsed);
        }
        BiasedBarcodeCollection result = new BiasedBarcodeCollection(errorBarcodesWithPositions, numCellsFilteredLowUMIs);
        return result;
	}

	/**
	 * Gets a reference to the UMI strings from the cache.  Has the side effect of populating the cache with additional
	 * strings.  This reduces total memory footprint by returning references to repeated strings instead of
	 * holding new objects for the same UMI over and over.
	 * @param umis A list of strings to get references to
	 * @param umiStringCache The cache of strings holding references.
	 */
	private Collection<String> getUMIsFromCache (final Collection<String> umis, final StringInterner umiStringCache) {
		List<String> result = new ArrayList<>(umis.size());
		for (String umi: umis)
			result.add(umiStringCache.intern(umi));
		return (result);
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
	public String fixUMI (final String cellBarcode, final String umi, final int errorPosition) {
		// 0 based, from end of cell barcode.
		int badBasesUMI=umi.length()-errorPosition;
		int lastBase = cellBarcode.length();
		int firstBaseToPad = lastBase-badBasesUMI-1;
		String cellBCBases=cellBarcode.substring(firstBaseToPad, cellBarcode.length());

		String umiRemaining=umi.substring(0, errorPosition-1);
		return cellBCBases+umiRemaining;
	}

	/**
	 * Set up the UMI Iterator.
	 * @return
	 */
	public UMIIterator prepareUMIIterator() {
		List<String> barcodes=getCellBarcodes();
		return new UMIIterator(SamFileMergeUtil.mergeInputs(INPUT, false, samReaderFactory),
                this.GENE_EXON_TAG, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, this.STRAND_TAG, this.READ_MQ,
                false, false, barcodes, true);
	}

	public List<String> getCellBarcodes () {

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

	BeadSynthesisErrorTypes getEnhancedErrorType (final BeadSynthesisErrorData data, final double extremeBaseRatio, final DetectPrimerInUMI detectPrimerTool) {
		BeadSynthesisErrorTypes errorType = data.getErrorType(extremeBaseRatio);
		//base case, error is not a single UMI.
		if (errorType!=BeadSynthesisErrorTypes.SINGLE_UMI)
			return errorType;
		// if there's a primer, run detection.
		if (detectPrimerTool!=null) {
			// a single UMI-style error, does the most common UMI match the primer?
			String singleUMI = data.getUMICounts().getKeysOrderedByCount(true).get(0);
			boolean primerDetected = detectPrimerTool.isStringInPrimer(singleUMI, this.EDIT_DISTANCE);
			if (primerDetected)
				return BeadSynthesisErrorTypes.PRIMER;
			else
				return errorType;
		}
		return errorType;
	}

	/**
	 * Picks a number of bases to pad.
	 * If errorPosition =-1, then don't pad any bases.
	 */
	String padCellBarcode (final String cellBarcode, final int errorPosition, final int umiLength) {
		if (errorPosition==-1) return (cellBarcode);

		// 0 based, from end of cell barcode.
		int badBasesUMI=umiLength-errorPosition;
		int lastBase = cellBarcode.length();
		int firstBaseToPad = lastBase-badBasesUMI-1;

		char [] charAr = cellBarcode.toCharArray();
		for (int i=firstBaseToPad; i<lastBase; i++)
			charAr[i]=this.PAD_CHARACTER;
		return new String (charAr);
	}


}
