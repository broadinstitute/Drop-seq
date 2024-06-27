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
package org.broadinstitute.dropseqrna.barnyard;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintWriter;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import picard.annotation.LocusFunction;


public class DigitalExpressionTest {


	// generate an input file with 5 barcodes and 3 genes.
	private static final File IN_CELL_BARCODE_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene.cellbarcodes.txt");
	private static final File IN_CELL_BARCODE_WITH_EXTRAS_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene_with_extras.cellbarcodes.txt");

	// expected results for standard coding strand specific DGE
	private static final File EXPECTED_OUTFILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene.dge.txt");
	private static final File EXPECTED_OUTFILE_SUMMARY = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene.dge_summary.txt");
	public static final File EXPECTED_OUTFILE_LONG = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene.dge_long.txt");

	private static final File EXPECTED_OUTPUT_SINGLE_CELL=new File("testdata/org/broadinstitute/transcriptome/barnyard/1_cell.dge.txt");
	//

	public static final String GENE_NAME_TAG="gn";
	public static final String GENE_STRAND_TAG="gs";
	public static final String GENE_FUNCTION_TAG="gf";
	public static final StrandStrategy STRAND_STRATEGY = StrandStrategy.SENSE;
	public static final List<LocusFunction> LOCUS_FUNCTION_LIST=new ArrayList<>(Arrays.asList(LocusFunction.CODING, LocusFunction.UTR));
	public static final String CELL_BARCODE_TAG="ZC";
	public static final String MOLECULAR_BARCODE_TAG = "XM";
	public static final int READ_MQ=10;
	private static final Boolean USE_STRAND_INFO=true;


	private static UMIIterator getUMIIterator (final File inFile) {
		List<String> cellBarcodes = Arrays.asList(DigitalExpressionTestUtil.barcodes);

		UMIIterator umiIterator = new UMIIterator(SamFileMergeUtil.mergeInputs(Collections.singletonList(inFile), false), GENE_NAME_TAG,
				GENE_STRAND_TAG, GENE_FUNCTION_TAG, STRAND_STRATEGY, LOCUS_FUNCTION_LIST, GeneFunctionCommandLineBase.DEFAULT_FUNCTIONAL_STRATEGY,
				CELL_BARCODE_TAG, MOLECULAR_BARCODE_TAG, READ_MQ, true, cellBarcodes, false, false,
				false, null);

		return (umiIterator);
	}

	// DigitalExpression I=5cell3gene_retagged.bam SUMMARY=5cell3gene.dge_summary.txt O=5cell3gene.dge.txt OUTPUT_LONG_FORMAT=5cell3gene.dge_long.txt CELL_BC_FILE=5cell3gene.cellbarcodes.txt

	@Test(dataProvider = "testDoWorkDataProvider")
	public void testDoWork (final File cellBCFile, final boolean omitMissingCells, final boolean testExpectedOutput) throws IOException {
		final File outFile = File.createTempFile("testDigitalExpression.", ".digital_expression.txt");
		final File summaryFile = File.createTempFile("testDigitalExpression.", ".digital_expression_summary.txt");
		final File cellBarcodesFile = File.createTempFile("testDigitalExpression.", ".selectedCellBarcodes.txt");
		final File longOutput = File.createTempFile("testDigitalExpression.", ".digital_expression_long.txt");
		outFile.deleteOnExit();
		summaryFile.deleteOnExit();
		cellBarcodesFile.deleteOnExit();
		longOutput.deleteOnExit();

		final DigitalExpression de = new DigitalExpression();
		de.INPUT = DigitalExpressionTestUtil.IN_FILE;
		de.CELL_BC_FILE = cellBCFile;
		de.OUTPUT = outFile;
		de.SUMMARY = summaryFile;
		de.OUTPUT_LONG_FORMAT=longOutput;
		de.OMIT_MISSING_CELLS = omitMissingCells;
		// the headers aren't going to match up because they contain specific path info.
		// de.UNIQUE_EXPERIMENT_ID = "test";

		int result = de.doWork();
		Assert.assertEquals(result, 0);

		if (testExpectedOutput) {
			Assert.assertTrue(FileUtils.contentEquals(outFile, EXPECTED_OUTFILE));
			Assert.assertTrue(FileUtils.contentEquals(summaryFile, EXPECTED_OUTFILE_SUMMARY));
			Assert.assertTrue(FileUtils.contentEquals(longOutput, EXPECTED_OUTFILE_LONG));
		}
	}

	@DataProvider(name = "testDoWorkDataProvider")
	public Object[][] testDoWorkDataProvider() {
		return new Object[][] {
				{IN_CELL_BARCODE_FILE, false, true},
				{IN_CELL_BARCODE_FILE, true, true},
				{IN_CELL_BARCODE_WITH_EXTRAS_FILE, true, true},
				{IN_CELL_BARCODE_WITH_EXTRAS_FILE, false, false}
		};
	}

	//TODO: set up the proper output files.
	@Test (enabled=true)
	public void testDoWorkSingleBarcode () {
		File outFile=null;
		File cellBarcodesFile=null;
		try {
			outFile = File.createTempFile("testDigitalExpression.", ".digital_expression.txt");
	        cellBarcodesFile = File.createTempFile("testDigitalExpression.", ".selectedCellBarcodes.txt");
	        outFile.deleteOnExit();
	        cellBarcodesFile.deleteOnExit();
		} catch (IOException e) {
			e.printStackTrace();
		}

        final DigitalExpression de = new DigitalExpression();
		de.INPUT = DigitalExpressionTestUtil.IN_FILE;
		de.NUM_CORE_BARCODES=1;
		de.OUTPUT = outFile;
        // the headers aren't going to match up because they contain specific path info.
        // de.UNIQUE_EXPERIMENT_ID = "test";

        int result = de.doWork();
        Assert.assertEquals(result, 0);

        try {
			Assert.assertTrue (FileUtils.contentEquals(outFile, EXPECTED_OUTPUT_SINGLE_CELL));
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	@Test ()
	public void testCustomCommandLineValidation1 () {
		// if header is set, then must have a UEI.
		final DigitalExpression de = new DigitalExpression();
		de.OUTPUT_HEADER=true;
        de.UNIQUE_EXPERIMENT_ID="foo";
        String [] o = de.customCommandLineValidation();
        Assert.assertNull(o);

        de.UNIQUE_EXPERIMENT_ID=null;
        o = de.customCommandLineValidation();
        Assert.assertNotNull(o);

        de.OUTPUT_HEADER=null;
        o = de.customCommandLineValidation();
        Assert.assertNull(o);



	}

	@Test(groups={"dropseq", "transcriptome"})
	public void DGEIntegrationTest () {
		UMIIterator u = getUMIIterator (DigitalExpressionTestUtil.IN_FILE);
		UMICollection batch;
		int count=0;

		while ((batch=u.next())!=null) {
			if (batch==null || batch.isEmpty())
				continue;
			count++;
			String currentGene = batch.getGeneName();
			String currentCell = batch.getCellBarcode();
			int expectedReadcount = getReadCountsNew(currentGene, currentCell);
			int dgeReadCount = batch.getDigitalExpression(1, 0, true);
			Assert.assertEquals(dgeReadCount,expectedReadcount);

			int dgeNoCollapseExpected=getDGEWithoutCollapseNew(currentGene, currentCell);
			int dgeNoCollapseActual = batch.getDigitalExpression(1, 0, false);
			Assert.assertEquals(dgeNoCollapseActual,dgeNoCollapseExpected);

		}
		// assert that you actually looked at a batch of data!
		Assert.assertNotEquals(count,0);
	}

	public void testTwoGenesOnSameStrand () {
		List barcodes = Collections.singletonList("FOO");
		File inFile = new File ("");

		UMIIterator umiIterator = new UMIIterator(SamFileMergeUtil.mergeInputs(Collections.singletonList(inFile), false), GENE_NAME_TAG,
				GENE_STRAND_TAG, GENE_FUNCTION_TAG, this.STRAND_STRATEGY, this.LOCUS_FUNCTION_LIST, GeneFunctionCommandLineBase.DEFAULT_FUNCTIONAL_STRATEGY,
				this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, this.READ_MQ, true, barcodes);
	}





	//GET COUNTS OF READS
	// ~/samtools view -q 10 5cell3gene.bam |grep ZC:Z:ATCAGGGACAGA |grep HUMAN_3:42642106-42690227:NKTR  | sed -E 's/.*XM
	private int getReadCounts (final String gene, final String cell) {

		if (cell.equals("ATCAGGGACAGA") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 845;
		if (cell.equals("ATCAGGGACAGA") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 200;
		if (cell.equals("ATCAGGGACAGA") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 225;

		if (cell.equals("AGGGAAAATTGA") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 428;
		if (cell.equals("AGGGAAAATTGA") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 212;
		if (cell.equals("AGGGAAAATTGA") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 238;

		if (cell.equals("TTGCCTTACGCG") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 473;
		if (cell.equals("TTGCCTTACGCG") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 62;
		if (cell.equals("TTGCCTTACGCG") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 581;

		if (cell.equals("TGGCGAAGAGAT") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 612;
		if (cell.equals("TGGCGAAGAGAT") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 12;
		if (cell.equals("TGGCGAAGAGAT") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 166;

		if (cell.equals("TACAATTAAGGC") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 385;
		if (cell.equals("TACAATTAAGGC") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 160;
		if (cell.equals("TACAATTAAGGC") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 294;
		return -1;

	}

	//GET COUNTS OF READS
		// ~/samtools view -q 10 5cell3gene.bam |grep ZC:Z:ATCAGGGACAGA |grep HUMAN_3:42642106-42690227:NKTR  | sed -E 's/.*XM
		private int getReadCountsNew (final String gene, final String cell) {

			if (cell.equals("ATCAGGGACAGA") && gene.equals("NKTR") ) return 845;
			if (cell.equals("ATCAGGGACAGA") && gene.equals("CHUK") ) return 200;
			if (cell.equals("ATCAGGGACAGA") && gene.equals("SNRPA1") ) return 225;

			if (cell.equals("AGGGAAAATTGA") && gene.equals("NKTR") ) return 428;
			if (cell.equals("AGGGAAAATTGA") && gene.equals("CHUK") ) return 212;
			if (cell.equals("AGGGAAAATTGA") && gene.equals("SNRPA1") ) return 238;

			if (cell.equals("TTGCCTTACGCG") && gene.equals("NKTR") ) return 473;
			if (cell.equals("TTGCCTTACGCG") && gene.equals("CHUK") ) return 62;
			if (cell.equals("TTGCCTTACGCG") && gene.equals("SNRPA1") ) return 581;

			if (cell.equals("TGGCGAAGAGAT") && gene.equals("NKTR") ) return 612;
			if (cell.equals("TGGCGAAGAGAT") && gene.equals("CHUK") ) return 12;
			if (cell.equals("TGGCGAAGAGAT") && gene.equals("SNRPA1") ) return 166;

			if (cell.equals("TACAATTAAGGC") && gene.equals("NKTR") ) return 385;
			if (cell.equals("TACAATTAAGGC") && gene.equals("CHUK") ) return 160;
			if (cell.equals("TACAATTAAGGC") && gene.equals("SNRPA1") ) return 294;
			return -1;

		}

	// GET DGE (WITHOUT BC COLLAPSE)
	// ~/samtools view -q 10 5cell3gene.bam |grep ZC:Z:ATCAGGGACAGA |grep HUMAN_3:42642106-42690227:NKTR  | sed -E 's/.*XM:Z:([ACGT]+).*/\1/' |sort |uniq -c |wc -l
	private int getDGEWithoutCollapse (final String gene, final String cell) {

		if (cell.equals("ATCAGGGACAGA") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 84;
		if (cell.equals("ATCAGGGACAGA") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 19;
		if (cell.equals("ATCAGGGACAGA") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 34;

		if (cell.equals("AGGGAAAATTGA") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 54;
		if (cell.equals("AGGGAAAATTGA") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 14;
		if (cell.equals("AGGGAAAATTGA") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 45;

		if (cell.equals("TTGCCTTACGCG") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 40;
		if (cell.equals("TTGCCTTACGCG") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 4;
		if (cell.equals("TTGCCTTACGCG") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 92;

		if (cell.equals("TGGCGAAGAGAT") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 59;
		if (cell.equals("TGGCGAAGAGAT") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 2;
		if (cell.equals("TGGCGAAGAGAT") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 32;

		if (cell.equals("TACAATTAAGGC") && gene.equals("HUMAN_3:42642106-42690227:NKTR") ) return 40;
		if (cell.equals("TACAATTAAGGC") && gene.equals("HUMAN_10:101948055-101989376:CHUK") ) return 11;
		if (cell.equals("TACAATTAAGGC") && gene.equals("HUMAN_15:101821715-101835487:SNRPA1") ) return 32;
		return -1;

	}

	// GET DGE (WITHOUT BC COLLAPSE)
		// ~/samtools view -q 10 5cell3gene.bam |grep ZC:Z:ATCAGGGACAGA |grep HUMAN_3:42642106-42690227:NKTR  | sed -E 's/.*XM:Z:([ACGT]+).*/\1/' |sort |uniq -c |wc -l
		private int getDGEWithoutCollapseNew (final String gene, final String cell) {

			if (cell.equals("ATCAGGGACAGA") && gene.equals("NKTR") ) return 84;
			if (cell.equals("ATCAGGGACAGA") && gene.equals("CHUK") ) return 19;
			if (cell.equals("ATCAGGGACAGA") && gene.equals("SNRPA1") ) return 34;

			if (cell.equals("AGGGAAAATTGA") && gene.equals("NKTR") ) return 54;
			if (cell.equals("AGGGAAAATTGA") && gene.equals("CHUK") ) return 14;
			if (cell.equals("AGGGAAAATTGA") && gene.equals("SNRPA1") ) return 45;

			if (cell.equals("TTGCCTTACGCG") && gene.equals("NKTR") ) return 40;
			if (cell.equals("TTGCCTTACGCG") && gene.equals("CHUK") ) return 4;
			if (cell.equals("TTGCCTTACGCG") && gene.equals("SNRPA1") ) return 92;

			if (cell.equals("TGGCGAAGAGAT") && gene.equals("NKTR") ) return 59;
			if (cell.equals("TGGCGAAGAGAT") && gene.equals("CHUK") ) return 2;
			if (cell.equals("TGGCGAAGAGAT") && gene.equals("SNRPA1") ) return 32;

			if (cell.equals("TACAATTAAGGC") && gene.equals("NKTR") ) return 40;
			if (cell.equals("TACAATTAAGGC") && gene.equals("CHUK") ) return 11;
			if (cell.equals("TACAATTAAGGC") && gene.equals("SNRPA1") ) return 32;
			return -1;

		}


	@Test
	public void testDigitalExpression() throws IOException {
		DigitalExpressionTestUtil.makeDigitalExpressionFile(true);
	}

    /**
     * Have to go through CLP in order to make sure this works, because UEI check was in CLP validation.
     */
    @Test
    //TODO what is this test supposed to do?
    public void testDigitalExpressionWithHeaderButNoUei() throws IOException {
        final File[] tempFiles = prepareToDigitalExpression();
        final String[] dgeArgs = new String[]{
                "INPUT=" + DigitalExpressionTestUtil.IN_FILE.getAbsolutePath(),
                "OUTPUT=" + tempFiles[1].getAbsolutePath(),
                "SUMMARY=" + tempFiles[2].getAbsolutePath(),
                "CELL_BC_FILE=" + tempFiles[0].getAbsolutePath(),
                "OUTPUT_HEADER=false"
        };
        Assert.assertEquals(new DigitalExpression().instanceMain(dgeArgs), 0);
    }


	/**
     * Creates barcodes file, DGE outfile (empty) and DGE summary file (empty), in temp directory and marked as
     * deleteOnExit().
     * @return [cellBarcodeFile, DgeOutFile, DgeSummaryFile]
     * @throws IOException
     */
    private static File[] prepareToDigitalExpression() throws IOException {
        final File outFile = File.createTempFile("testDigitalExpression.", ".digital_expression.txt.gz");
        final File summaryFile = File.createTempFile("testDigitalExpression.", ".digital_expression_summary.txt");
        final File cellBarcodesFile = File.createTempFile("testDigitalExpression.", ".selectedCellBarcodes.txt");
        outFile.deleteOnExit();
        summaryFile.deleteOnExit();
        cellBarcodesFile.deleteOnExit();
        final ErrorCheckingPrintWriter writer = new ErrorCheckingPrintWriter(cellBarcodesFile);
        for (final String cellBarcode : DigitalExpressionTestUtil.barcodes)
			writer.println(cellBarcode);
        writer.close();
        return new File[]{cellBarcodesFile, outFile, summaryFile};
    }

    private static File STRAND_AND_FUNCTION_TESTDATA = new File("testdata/org/broadinstitute/dropseq/barnyard/DgeStrandFuncTest");

	//TODO: I have no idea what this test is trying to accomplish.  The strategies must be set to some non-null value.
	@Test(dataProvider = "testSuppressingStrandAndFunctionDataProvider", enabled = false)
	public void testSuppressingStrandAndFunction(final boolean suppressStrandStrategy,
												 final boolean suppressFunction,
												 final boolean fakeTagNames) throws IOException {
    	final String expectedResultLabel;
    	if (suppressFunction) {
    		if (suppressStrandStrategy) {
    			expectedResultLabel = "neither";
			} else {
    			expectedResultLabel = "strand";
			}
		} else if (suppressStrandStrategy) {
    		expectedResultLabel = "func";
		} else {
    		expectedResultLabel = "both";
		}
		final String summaryExtension = ".digital_expression_summary.txt";
		final File expectedResultFile = new File(STRAND_AND_FUNCTION_TESTDATA, expectedResultLabel + summaryExtension);
    	final DigitalExpression digitalExpression = new DigitalExpression();
    	digitalExpression.INPUT = new File(STRAND_AND_FUNCTION_TESTDATA, "DgeStrandFuncTest.bam");
    	digitalExpression.CELL_BC_FILE = new File(STRAND_AND_FUNCTION_TESTDATA, "DgeStrandFuncTest.cell_barcodes");
    	digitalExpression.SUMMARY = File.createTempFile("DgeStrandFuncTest.", summaryExtension);
    	digitalExpression.SUMMARY.deleteOnExit();
    	digitalExpression.OUTPUT = new File("/dev/null");
    	if (suppressStrandStrategy) {
    		digitalExpression.STRAND_STRATEGY = null;
    		if (fakeTagNames) {
    			digitalExpression.GENE_STRAND_TAG = "00";
			}
		}
    	if (suppressFunction) {
    		digitalExpression.LOCUS_FUNCTION_LIST.clear();
			if (fakeTagNames) {
				digitalExpression.GENE_FUNCTION_TAG = "00";
			}
		}
    	Assert.assertEquals(digitalExpression.doWork(), 0);
    	Assert.assertTrue(TestUtils.testMetricsFilesEqual(expectedResultFile, digitalExpression.SUMMARY));
	}

	@DataProvider(name = "testSuppressingStrandAndFunctionDataProvider")
	public Object[][] testSuppressingStrandAndFunctionDataProvider() {
    	final Object[][] ret = new Object[8][];
    	final boolean[] tf = {true, false};
    	int i = 0;
    	for (boolean a : tf) for (boolean b : tf) for (boolean c: tf) {
    		ret[i++] = new Object[]{a,b,c};
		}
    	return ret;
	}
}
