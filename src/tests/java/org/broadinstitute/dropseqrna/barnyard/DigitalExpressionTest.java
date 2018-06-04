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

import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintWriter;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.annotation.LocusFunction;

import java.io.File;
import java.io.IOException;
import java.util.*;


public class DigitalExpressionTest {


	// generate an input file with 5 barcodes and 3 genes.
	/**
	 * For the integration test, make a relatively small BAM file and read it in.
	 * Manually test the counts of the reads across UMIs - umi testing to collapse molecular barcodes is in it's own test.
	 *
	 * Data taken from 9-27-14 100 cells data set.
	 *
	 * cells
	 * ATCAGGGACAGA
	 * AGGGAAAATTGA
	 * TTGCCTTACGCG
	 * TGGCGAAGAGAT
	 * TACAATTAAGGC
	 *
	 * genes
	 * HUMAN_10:101948055-101989376:CHUK
	 * HUMAN_15:101821715-101835487:SNRPA1
	 * HUMAN_3:42642106-42690227:NKTR
	 *
	 * /fg/software/gap/gap_analysis/FilterBAMByTag I=100cells_star_bq10_noPseudoGenes.bam O=test.bam TAG=ZC TAG_VALUES_FILE=bc.txt ACCEPT_TAG=true
	 * samtools view -H test.bam > 5cell3gene.sam
	 * samtools view test.bam |grep -f genes.txt >> 5cell3gene.sam
	 * samtools view -Sb 5cell3gene.sam > 5cell3gene.bam
	 *
	 */

	private static final File IN_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene_retagged.bam");

	private String GENE_NAME_TAG="gn";
	private String GENE_STRAND_TAG="gs";
	private String GENE_FUNCTION_TAG="gf";
	private StrandStrategy STRAND_STRATEGY = StrandStrategy.SENSE;
	private List<LocusFunction> LOCUS_FUNCTION_LIST=new ArrayList<>(Arrays.asList(LocusFunction.CODING, LocusFunction.UTR));
	private String CELL_BARCODE_TAG="ZC";
	private String MOLECULAR_BARCODE_TAG = "XM";
	private int READ_MQ=10;
	private Boolean USE_STRAND_INFO=true;
    private static final String [] barcodes ={"ATCAGGGACAGA", "AGGGAAAATTGA", "TTGCCTTACGCG", "TGGCGAAGAGAT", "TACAATTAAGGC"};


	private UMIIterator getUMIIterator (final File inFile) {
		List<String> cellBarcodes = Arrays.asList(barcodes);

		UMIIterator umiIterator = new UMIIterator(SamFileMergeUtil.mergeInputs(Collections.singletonList(inFile), false), GENE_NAME_TAG,
				GENE_STRAND_TAG, GENE_FUNCTION_TAG, this.STRAND_STRATEGY, this.LOCUS_FUNCTION_LIST, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG,
        		this.READ_MQ, true, cellBarcodes);

		return (umiIterator);
	}


	@Test(groups={"dropseq", "transcriptome"})
	public void DGEIntegrationTest () {
		UMIIterator u = getUMIIterator (IN_FILE);
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
				GENE_STRAND_TAG, GENE_FUNCTION_TAG, this.STRAND_STRATEGY, this.LOCUS_FUNCTION_LIST, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG,
        		this.READ_MQ, true, barcodes);
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
		makeDigitalExpressionFile(true);
	}

    /**
     * Have to go through CLP in order to make sure this works, because UEI check was in CLP validation.
     */
    @Test
    //TODO what is this test supposed to do?
    public void testDigitalExpressionWithHeaderButNoUei() throws IOException {
        final File[] tempFiles = prepareToDigitalExpression();
        final String[] dgeArgs = new String[]{
                "INPUT=" + IN_FILE.getAbsolutePath(),
                "OUTPUT=" + tempFiles[1].getAbsolutePath(),
                "SUMMARY=" + tempFiles[2].getAbsolutePath(),
                "CELL_BC_FILE=" + tempFiles[0].getAbsolutePath(),
                "OUTPUT_HEADER=false"
        };
        Assert.assertEquals(new DigitalExpression().instanceMain(dgeArgs), 0);
    }

	/**
	 * @return A digital expression file, which is marked as deleteOnExit
     */
    public static File makeDigitalExpressionFile(final boolean outputHeader) throws IOException {
		final File outFile = File.createTempFile("testDigitalExpression.", ".digital_expression.txt.gz");
		final File summaryFile = File.createTempFile("testDigitalExpression.", ".digital_expression_summary.txt");
        final File cellBarcodesFile = File.createTempFile("testDigitalExpression.", ".selectedCellBarcodes.txt");
		outFile.deleteOnExit();
		summaryFile.deleteOnExit();
        cellBarcodesFile.deleteOnExit();
        final ErrorCheckingPrintWriter writer = new ErrorCheckingPrintWriter(cellBarcodesFile);
        for (final String cellBarcode : barcodes)
			writer.println(cellBarcode);
        writer.close();
        final DigitalExpression de = new DigitalExpression();
		de.INPUT = IN_FILE;
		de.OUTPUT = outFile;
		de.SUMMARY = summaryFile;
        de.CELL_BC_FILE = cellBarcodesFile;
        de.OUTPUT_HEADER = outputHeader;
        de.UNIQUE_EXPERIMENT_ID = "UIE" + new Random().nextInt();
        Assert.assertEquals(de.doWork(), 0);
        return outFile;
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
        for (final String cellBarcode : barcodes)
			writer.println(cellBarcode);
        writer.close();
        return new File[]{cellBarcodesFile, outFile, summaryFile};
    }

}
