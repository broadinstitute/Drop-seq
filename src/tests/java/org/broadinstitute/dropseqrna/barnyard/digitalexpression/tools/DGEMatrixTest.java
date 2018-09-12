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
package org.broadinstitute.dropseqrna.barnyard.digitalexpression.tools;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import org.junit.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.util.TestUtil;

public class DGEMatrixTest {

	private final File exampleOne = new File("testdata/org/broadinstitute/transcriptome/barnyard/digitalexpression/dge_example1.txt.gz");
	private final File exampleTwo = new File("testdata/org/broadinstitute/transcriptome/barnyard/digitalexpression/dge_example2.txt.gz");
	private final File exampleThree = new File("testdata/org/broadinstitute/transcriptome/barnyard/digitalexpression/dge_example3.txt.gz");
	private final File preMergedExample = new File("testdata/org/broadinstitute/transcriptome/barnyard/digitalexpression/dge_example_merged.txt.gz");
	private final File preMergedExample2 = new File("testdata/org/broadinstitute/transcriptome/barnyard/digitalexpression/dge_example_merged2.txt.gz");
	private final File matrixMarketMatrixFile = new File ("testdata/org/broadinstitute/transcriptome/barnyard/digitalexpression/tenXMatrixMarket.mtx");
	private final File matrixMarketCellFile = new File ("testdata/org/broadinstitute/transcriptome/barnyard/digitalexpression/tenXMatrixMarketCellBarcodes.tsv");
	private final File matrixMarketGeneFile = new File ("testdata/org/broadinstitute/transcriptome/barnyard/digitalexpression/tenXMatrixMarketGenes.tsv");



	@Test
	public void parseFile() {
		final DGEMatrix result= DGEMatrix.parseFile(exampleOne);

		// test cells the same set, in the same order.
		final List<String> cells = result.getCellBarcodes();
		final String [] expectedCells = {"CAATCCGACAAC", "CACTAAAGCCAG", "TCCCTTCAAGTA", "ATGGTCTCAAAC", "CCTTCCATGCGA"};
		final String [] actualCells = cells.toArray(new String [cells.size()]);
		Assert.assertArrayEquals(expectedCells, actualCells);

		// test getting genes out.  They are returned in alphabetical order.
		final String [] expectedGenes ={"Arl6ip1", "Prkg1", "Rasa1", "Exosc2", "Gimap5","Ranbp6","Tet3", "Nmnat3"};
		final List<String> genes = result.getGenes();
		final String [] actualGenes = genes.toArray(new String [genes.size()]);
		Assert.assertArrayEquals(expectedGenes, actualGenes);

		// test the expression data.  Test a few genes.
		double [] expression = result.getExpression("Arl6ip1");
		final double [] expected = {25,30,36,14,29};
		for (int i=0; i<expected.length; i++)
			Assert.assertEquals(expected[i], expression[i], 0.001);

		expression = result.getExpression("Rasa1");
		final double [] expected2 = {8,6,4,4,4};
		for (int i=0; i<expected.length; i++)
			Assert.assertEquals(expected2[i], expression[i], 0.001);

		expression = result.getExpression("Nmnat3");
		final double [] expected3 = {0,0,3,0,0};
		for (int i=0; i<expected.length; i++)
			Assert.assertEquals(expected3[i], expression[i], 0.001);
	}

	@Test
	public void testMerge () {
		final DGEMatrix resultOne= DGEMatrix.parseFile(exampleOne, "ExpA_");
		final DGEMatrix resultTwo= DGEMatrix.parseFile(exampleTwo, "ExpB_");

		final DGEMatrix result = resultOne.merge(resultTwo);
		final DGEMatrix preMerged= DGEMatrix.parseFile(preMergedExample);

		// compare cell barcodes
		final List<String> cells = result.getCellBarcodes();
		final List<String> newCells = preMerged.getCellBarcodes();
		Assert.assertArrayEquals(cells.toArray(new String [cells.size()]), newCells.toArray(new String [newCells.size()]));

		// compare genes
		final List<String> genes = result.getGenes();
		final List<String> newGenes = preMerged.getGenes();
		Assert.assertArrayEquals(genes.toArray(new String [genes.size()]), newGenes.toArray(new String [newGenes.size()]));

		// compare expression
		final double [] [] expectedExp = result.getExpressionMatrix();
		final double [] [] obsExp = result.getExpressionMatrix();

		//dimensions
		Assert.assertEquals(expectedExp.length, obsExp.length);
		Assert.assertEquals(expectedExp[0].length, obsExp[0].length);

		for (int row=0; row<expectedExp.length; row++)
			for (int col=0; col<expectedExp[0].length; col++)
				Assert.assertEquals(expectedExp[row][col], obsExp[row][col], 0.001);

	}

	@Test
	// merge a file with itself, so the result is double the input.
	public void testMergeWithCollapseSimple () {
		final DGEMatrix resultOne= DGEMatrix.parseFile(exampleOne);
		final DGEMatrix resultTwo= DGEMatrix.parseFile(exampleOne);

		final DGEMatrix result = resultOne.mergeWithCollapse(resultTwo);
		for (String gene: result.getGenes()) {
			double [] expression = result.getExpression(gene);
			double [] expressionExpected = resultTwo.getExpression(gene);
			Assert.assertEquals(expression.length, expressionExpected.length);
			for (int i=0; i<expression.length; i++) {
				double e = expression[i];
				double expected = expressionExpected[i]*2;
				Assert.assertEquals(e, expected, 0.1);
			}
		}
	}

	@Test
	// A harder merge where not all cell barcodes are duplicated.
	public void testMergeWithCollapseHarder () {
		final DGEMatrix resultOne= DGEMatrix.parseFile(exampleOne);
		final DGEMatrix resultTwo= DGEMatrix.parseFile(exampleThree);
		final DGEMatrix preMerged= DGEMatrix.parseFile(preMergedExample2);

		final DGEMatrix result = resultOne.mergeWithCollapse(resultTwo);

		Assert.assertEquals(result.getGenes().size(), preMerged.getGenes().size());
		Assert.assertEquals(result.getCellBarcodes().size(), preMerged.getCellBarcodes().size());

		for (String gene: result.getGenes()) {
			double [] expression = result.getExpression(gene);
			double [] expressionExpected = preMerged.getExpression(gene);
			Assert.assertEquals(expression.length, expressionExpected.length);
			for (int i=0; i<expression.length; i++) {
				double e = expression[i];
				double expected = expressionExpected[i];
				Assert.assertEquals(e, expected, 0.1);
			}
		}
	}

	@Test
	// A harder merge where there's no cell barcode overlap, so should be the same as the merge() result.
	public void testMergeWithCollapseHarder2 () {
		final DGEMatrix resultOne= DGEMatrix.parseFile(exampleOne);
		final DGEMatrix resultTwo= DGEMatrix.parseFile(exampleTwo);
		final DGEMatrix preMerged= DGEMatrix.parseFile(preMergedExample);

		final DGEMatrix result = resultOne.mergeWithCollapse(resultTwo);

		Assert.assertEquals(result.getGenes().size(), preMerged.getGenes().size());
		Assert.assertEquals(result.getCellBarcodes().size(), preMerged.getCellBarcodes().size());

		for (String gene: result.getGenes()) {
			double [] expression = result.getExpression(gene);
			double [] expressionExpected = preMerged.getExpression(gene);
			Assert.assertEquals(expression.length, expressionExpected.length);
			for (int i=0; i<expression.length; i++) {
				double e = expression[i];
				double expected = expressionExpected[i];
				Assert.assertEquals(e, expected, 0.1);
			}
		}
	}

	@Test
	public void testCellBarcodePrefix () {
		final DGEMatrix result= DGEMatrix.parseFile(exampleOne, "PREFIX_");
		final List<String> cells = result.getCellBarcodes();
		final String [] expectedCells = {"PREFIX_CAATCCGACAAC", "PREFIX_CACTAAAGCCAG", "PREFIX_TCCCTTCAAGTA", "PREFIX_ATGGTCTCAAAC", "PREFIX_CCTTCCATGCGA"};
		final String [] actualCells = cells.toArray(new String [cells.size()]);
		Assert.assertArrayEquals(expectedCells, actualCells);
	}

	@Test
	public void instantiateFrom2DMatrix() {
		final String [] expectedCells = {"CAATCCGACAAC", "CACTAAAGCCAG", "TCCCTTCAAGTA", "ATGGTCTCAAAC", "CCTTCCATGCGA"};
		final List<String> barcodes = Arrays.asList(expectedCells);
		final String [] expectedGenes ={"Arl6ip1", "Exosc2", "Gimap5", "Nmnat3", "Prkg1","Ranbp6","Rasa1", "Tet3"};
		final List<String> genesList=Arrays.asList(expectedGenes);

		final double [] [] data ={{25,30,36,14,29},{0,0,0,0,1},{5,0,0,0,0},{0,0,3,0,0},{1,0,0,0,0},{1,3,3,2,1},{8,6,4,4,4},{1,2,1,1,0}};
		final DGEMatrix result = new DGEMatrix(barcodes, genesList, data);

		// test cells the same set, in the same order.
		final List<String> cells = result.getCellBarcodes();
		final String [] actualCells = cells.toArray(new String [cells.size()]);
		Assert.assertArrayEquals(expectedCells, actualCells);

		// test getting genes out.  They are returned in alphabetical order.
		final List<String> genes = result.getGenes();
		final String [] actualGenes = genes.toArray(new String [genes.size()]);
		Assert.assertArrayEquals(expectedGenes, actualGenes);

		// test the expression data.  Test a few genes.
		double [] expression = result.getExpression("Arl6ip1");
		final double [] expected = {25,30,36,14,29};
		for (int i=0; i<expected.length; i++)
			Assert.assertEquals(expected[i], expression[i], 0.001);

		expression = result.getExpression("Rasa1");
		final double[] expected2 = { 8, 6, 4, 4, 4 };
		for (int i = 0; i < expected.length; i++)
			Assert.assertEquals(expected2[i], expression[i], 0.001);

		expression = result.getExpression("Nmnat3");
		final double[] expected3 = { 0, 0, 3, 0, 0 };
		for (int i = 0; i < expected.length; i++)
			Assert.assertEquals(expected3[i], expression[i], 0.001);

		final double [] [] exp = result.getExpressionMatrix();
		for (int row=0; row<exp.length; row++)
			for (int col=0; col<exp[0].length; col++)
				Assert.assertEquals(data[row][col], exp[row][col], 0.001);

		// test null expression
		final double [] expNull = result.getExpression("FAKE");
		Assert.assertNull(expNull);
	}

	@Test
	public void testRemoveCells() {
		final DGEMatrix result= DGEMatrix.parseFile(exampleOne);
		final String [] cellsToRemove = {"CACTAAAGCCAG", "ATGGTCTCAAAC"};
		result.removeCellBarcodes(Arrays.asList(cellsToRemove));

		// test expression
		double [] expression = result.getExpression("Rasa1");
		final double [] expected = { 8, 4, 4 };
		for (int i = 0; i < expected.length; i++)
			Assert.assertEquals(expected[i], expression[i], 0.001);

		expression = result.getExpression("Exosc2");
		final double [] expected2 = {0,0,1};
		for (int i = 0; i < expected2.length; i++)
			Assert.assertEquals(expected2[i], expression[i], 0.001);

		expression = result.getExpression("Tet3");
		final double [] expected3 = {1,1,0};
		for (int i = 0; i < expected3.length; i++)
			Assert.assertEquals(expected3[i], expression[i], 0.001);

	}

	@Test
	public void testRemoveGenes() {
		final DGEMatrix result= DGEMatrix.parseFile(exampleOne);

		// test the expression data.
		final double [] expression = result.getExpression("Arl6ip1");
		final double [] expected = {25,30,36,14,29};
		for (int i=0; i<expected.length; i++)
			Assert.assertEquals(expected[i], expression[i], 0.001);
		Assert.assertTrue(result.getGenes().contains("Arl6ip1"));

		final String [] genesToRemove = {"Arl6ip1"};
		result.removeGenes(Arrays.asList(genesToRemove));

		Assert.assertNull(result.getExpression("Arl6ip1"));
		Assert.assertFalse(result.getGenes().contains("Arl6ip1"));
	}

	@Test
	public void testRemoveGenes2() {
		final DGEMatrix result= DGEMatrix.parseFile(exampleOne);

		// test the expression data.
		final String [] genesToRemove = {"Nmnat3", "Prkg1", "Gimap5"};
		result.removeGenes(Arrays.asList(genesToRemove));

		Assert.assertNull(result.getExpression("Nmnat3"));
		Assert.assertNull(result.getExpression("Prkg1"));
		Assert.assertNull(result.getExpression("Gimap5"));
		Assert.assertTrue(result.getGenes().contains("Arl6ip1"));
		Assert.assertFalse(result.getGenes().contains("Nmnat3"));
		Assert.assertFalse(result.getGenes().contains("Prkg1"));
	}



	@Test
	public void testRemoveCellsWithLowExpression() {
		/*
		 *  apply (a, 2, function (x) length(which(x>0)))
		 *	CAATCCGACAAC CACTAAAGCCAG TCCCTTCAAGTA ATGGTCTCAAAC CCTTCCATGCGA
         * 	 		6            4            5            4            4
		 */
		final DGEMatrix result= DGEMatrix.parseFile(exampleOne);
		final int numGenes=4;
		result.removeCellsWithLowExpression(numGenes);
		List<String> cells = result.getCellBarcodes();
		final String [] actualCells = cells.toArray(new String [cells.size()]);
		String [] expectedCells ={"CAATCCGACAAC", "TCCCTTCAAGTA"};
		Assert.assertArrayEquals(expectedCells, actualCells);

		// test 2.
		final DGEMatrix result2= DGEMatrix.parseFile(exampleOne);
		final int numGenes2=3;
		result.removeCellsWithLowExpression(numGenes2);
		List<String> cells2 = result2.getCellBarcodes();
		final String [] actualCells2 = cells2.toArray(new String [cells2.size()]);
		String [] expectedCells2 ={"CAATCCGACAAC", "CACTAAAGCCAG", "TCCCTTCAAGTA", "ATGGTCTCAAAC","CCTTCCATGCGA"};
		Assert.assertArrayEquals(expectedCells, actualCells);

	}

	@Test
	public void testRemoveGenesWithLowExpression() {
		/*
		 * apply (a, 1, function (x) length(which(x>0)))
		 * Arl6ip1   Prkg1   Rasa1  Exosc2  Gimap5  Ranbp6    Tet3  Nmnat3
      			5       1       5       1       1       5       4       1
		 */
		final DGEMatrix result= DGEMatrix.parseFile(exampleOne);
		final int numCells=1;
		result.removeGenesWithLowExpression(numCells);
		List<String> genes = result.getGenes();
		final String [] actualGenes = genes.toArray(new String [genes.size()]);
		String [] expectedGenes ={"Arl6ip1", "Rasa1", "Ranbp6", "Tet3"};
		Assert.assertArrayEquals(expectedGenes, actualGenes);

	}

	@Test
	public void testMatrixMarketRoundTrip() throws IOException {
        final DGEMatrix original = DGEMatrix.parseFile(preMergedExample);
        final File mmFile = writeMatrixMarket(original);
        final DGEMatrix mm = DGEMatrix.parseFile(mmFile);
        assertDgesEqual(mm, original);
	}

    @Test
    public void testMatrixMarketWithPrefix() throws IOException {
        final DGEMatrix original = DGEMatrix.parseFile(preMergedExample);
        final File mmFile = writeMatrixMarket(original);
        final String prefix = "prefix_";
        final DGEMatrix mm = DGEMatrix.parseFile(mmFile, prefix);
        final DGEMatrix originalWithPrefix = DGEMatrix.parseFile(preMergedExample, prefix);
        assertDgesEqual(mm, originalWithPrefix);
    }

    @Test
    public void parse10xGenomicsMatrixMarketFile () {

    	final DGEMatrix d = DGEMatrix.parse10XGenomicsMatrixMarket(matrixMarketMatrixFile, matrixMarketGeneFile, matrixMarketCellFile, null);
    	d.toDenseMatrix();
    	double [] [] values =d.getExpressionMatrix();
    	double [] [] expected = {{1,0,0,4,0},{0,4,0,0,0}, {0,0,9,0,0}, {0,8,0,16,20},{0,0,0,0,25}};

    	for (int i=0; i<expected.length; i++)
			for (int j=0; j<expected[i].length; j++)
				Assert.assertEquals(expected[i][j], values[i][j], 0.001);

    }

    /*
    @Test
    public void testWriteDenseDGEFile () {
    	// round trip.

    	File temp=null;
    	try {
			temp = File.createTempFile("DenseDGE", ".txt.gz");
			temp.deleteOnExit();
		} catch (IOException e) {
			e.printStackTrace();
		}

    	DGEMatrix result= DGEMatrix.parseFile(exampleOne);
    	result.writeDenseDgeFile(temp, false);
    	DGEMatrix resultNew= DGEMatrix.parseFile(temp);
    	Assert.assertTrue(result.getCellBarcodes().equals(resultNew.getCellBarcodes()));
    	Assert.assertTrue(result.getGenes().equals(resultNew.getGenes()));

    	Assert.assertTrue(result.getExpressionMatrix().equals(resultNew.getExpressionMatrix()));



    }
	*/


    private File writeMatrixMarket(final DGEMatrix mat) throws IOException {
        final File tempDir = TestUtil.getTempDirectory("DGEMatrixTest.", ".tmp");
        tempDir.deleteOnExit();
        final File mmFile = File.createTempFile("DGEMatrixTest.", ".mm.gz", tempDir);
        mmFile.deleteOnExit();
        mat.writeFile(mmFile, true, DGEMatrix.FileFormat.MM_SPARSE);
        return mmFile;
    }

    private void assertDgesEqual(final DGEMatrix dge1, final DGEMatrix dge2) {
        Assert.assertEquals(dge1.getGenes(), dge2.getGenes());
        Assert.assertEquals(dge1.getCellBarcodes(), dge2.getCellBarcodes());
        Assert.assertEquals(dge1.getMatrix(), dge2.getMatrix());
    }


}
