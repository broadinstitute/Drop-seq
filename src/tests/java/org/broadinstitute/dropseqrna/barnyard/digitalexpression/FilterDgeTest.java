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
package org.broadinstitute.dropseqrna.barnyard.digitalexpression;

import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.dropseqrna.barnyard.DigitalExpression;
import org.broadinstitute.dropseqrna.barnyard.DigitalExpressionTest;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintWriter;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class FilterDgeTest {

	// cells CAATCCGACAAC	CACTAAAGCCAG	TCCCTTCAAGTA	ATGGTCTCAAAC	CCTTCCATGCGA
	// genes Arl6ip1 Prkg1 Rasa1 Exosc2 Gimap5 Ranbp6 Tet3 Nmnat3

	/*
	private final File exampleOne = new File("testdata/org/broadinstitute/dropseq/private/barnyard/digitalexpression/dge_example1.txt.gz");
	private final File retainCells=new File("testdata/org/broadinstitute/dropseq/private/barnyard/digitalexpression/retainCells.txt");
	private final File removeCells=new File("testdata/org/broadinstitute/dropseq/private/barnyard/digitalexpression/removeCells.txt");
	private final File retainGenes=new File("testdata/org/broadinstitute/dropseq/private/barnyard/digitalexpression/retainGenes.txt");
	private final File removeGenes=new File("testdata/org/broadinstitute/dropseq/private/barnyard/digitalexpression/removeGenes.txt");
	 */



    private static final String[] CELLS_TO_RETAIN = {"ATCAGGGACAGA", "TTGCCTTACGCG", "TACAATTAAGGC"};

	@Test(dataProvider = "testFilteringWithHeadersDataProvider")
	public void testFilteringWithHeaders(final boolean inputHasHeader, final boolean outputHasHeader) throws IOException {
		final File inDge = DigitalExpressionTest.makeDigitalExpressionFile(inputHasHeader);
        final File outDge = File.createTempFile("testFilteringWithHeaders.", DigitalExpressionTest.DIGITAL_EXPRESSION_EXTENSION);
        outDge.deleteOnExit();
        final File cellsFile = makeCellsToRetainFile();
        final FilterDge fd = new FilterDge();
        fd.INPUT = inDge;
        fd.OUTPUT = outDge;
        fd.OUTPUT_HEADER = outputHasHeader;
        fd.CELLS_RETAIN = cellsFile;
        Assert.assertEquals(fd.doWork(), 0);
	}

    private File makeCellsToRetainFile() throws IOException {
        final File cellsFile = File.createTempFile("testFilteringWithHeaders.", ".cells.list");
        cellsFile.deleteOnExit();
        final PrintWriter writer = new ErrorCheckingPrintWriter(cellsFile);
        for (final String cell : CELLS_TO_RETAIN)
			writer.println(cell);
        writer.close();
        return cellsFile;
    }

    @DataProvider(name = "testFilteringWithHeadersDataProvider")
    public Object[][] testFilteringWithHeadersDataProvider() {
        return new Object[][]{
                {false, false},
                {false, true},
                {true, false},
                {true, true}
        };
    }

    @Test(dataProvider = "testFilteringWithSummaryDataProvider")
    public void testFilteringWithSummary(final boolean useInputSummary) throws IOException {
	    final boolean inputHasHeader = false;
        final File inDge = DigitalExpressionTest.makeDigitalExpressionFile(inputHasHeader);
        final File inDgeSummary = DigitalExpressionTest.makeSummaryPathFromDgePath(inDge);
        final File outDge = File.createTempFile("testFilteringWithSummary.", DigitalExpressionTest.DIGITAL_EXPRESSION_EXTENSION);
        outDge.deleteOnExit();
        final File cellsFile = makeCellsToRetainFile();
        final FilterDge fd = new FilterDge();
        fd.INPUT = inDge;
        fd.OUTPUT = outDge;
        fd.OUTPUT_HEADER = false;
        fd.CELLS_RETAIN = cellsFile;
        fd.INPUT_SUMMARY = useInputSummary? inDgeSummary: null;
        fd.OUTPUT_SUMMARY = File.createTempFile("testFilteringWithSummary.", DigitalExpressionTest.DIGITAL_EXPRESSION_SUMMARY_EXTENSION);
        fd.OUTPUT_SUMMARY.deleteOnExit();
        Assert.assertEquals(fd.doWork(), 0);
        List<DigitalExpression.DESummary> summaries = MetricsFile.readBeans(fd.OUTPUT_SUMMARY);
        Assert.assertEquals(CELLS_TO_RETAIN.length, summaries.size());
        final Set<String> cellsToRetain = new HashSet<>(Arrays.asList(CELLS_TO_RETAIN));
        for (final DigitalExpression.DESummary summary : summaries) {
            Assert.assertTrue(cellsToRetain.contains(summary.CELL_BARCODE));
        }
    }

    @DataProvider(name = "testFilteringWithSummaryDataProvider")
    public Object[][] testFilteringWithSummaryDataProvider() {
	    return new Object[][] {
                {true},
                {false}
        };
    }
}
