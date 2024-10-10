/*
 * MIT License
 *
 * Copyright 2021 Broad Institute
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

import org.broadinstitute.dropseqrna.barnyard.MarkChimericReads;
import org.broadinstitute.dropseqrna.barnyard.ChimericUmi.CHIMERIC_STRATEGY;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.nio.PicardHtsPath;

import java.io.File;
import java.util.Arrays;

public class MarkChimericReadsTest {
	public static File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/barnyard");
	public static final File CELL_BC_FILE = new File(TEST_DATA_DIR, "selectedCellBarcodes.txt");
	public static final File EXPECTED_CHIMERIC_REPORT = new File(TEST_DATA_DIR, "umiReuse.chimeric_report.txt");

	@Test(enabled = true)
	public void testUmiReuse() {
		final MarkChimericReads clp = new MarkChimericReads();
		clp.OUTPUT = TestUtils.getTempReportFile("MarkChimericReads.", ".sam");
		clp.OUTPUT.deleteOnExit();
		clp.OUTPUT_REPORT = TestUtils.getTempReportFile("chimeric_report.", ".txt");
		clp.OUTPUT_REPORT.deleteOnExit();
		clp.METRICS = TestUtils.getTempReportFile("MarkChimericReads.", ".chimeric_read_metrics");
		clp.METRICS.deleteOnExit();
		clp.INPUT = Arrays.asList(new PicardHtsPath(new File(TEST_DATA_DIR, "MarkChimericReads.input.sam")));
		clp.CELL_BC_FILE = CELL_BC_FILE;
		Assert.assertEquals(clp.doWork(), 0);
		TestUtils.assertSamFilesSame(clp.OUTPUT, new File(TEST_DATA_DIR, "MarkChimericReads.umiReuse.sam"), false);
		Assert.assertTrue(TestUtils.testFilesSame(clp.OUTPUT_REPORT, EXPECTED_CHIMERIC_REPORT));
		Assert.assertTrue(TestUtils.testFilesSame(clp.METRICS, new File(TEST_DATA_DIR, "umiReuse.chimeric_read_metrics")));
	}

	// MarkChimericReads.input2.sam
	@Test(enabled = true)
	public void testUmiReuse10xStrategy() {
		// in this case, for molecular barcode ATTTAGTAAATG, there are 2 reads for USMG5 but only one for CXorf56.
		// USMG5 is retained, CXorf56 is marked as chimeric.
		// This is the 10x strategy for removing chimeric reads
		final MarkChimericReads clp = new MarkChimericReads();
		clp.OUTPUT = TestUtils.getTempReportFile("MarkChimericReads.", ".sam");
		clp.OUTPUT.deleteOnExit();
		clp.OUTPUT_REPORT = TestUtils.getTempReportFile("chimeric_report.", ".txt");
		clp.OUTPUT_REPORT.deleteOnExit();
		clp.METRICS = TestUtils.getTempReportFile("MarkChimericReads.", ".chimeric_read_metrics");
		clp.METRICS.deleteOnExit();
		clp.INPUT = Arrays.asList(new PicardHtsPath(new File(TEST_DATA_DIR, "MarkChimericReads.input2.sam")));
		clp.CELL_BC_FILE = new File(TEST_DATA_DIR, "selectedCellBarcodes.txt");
		clp.STRATEGY = CHIMERIC_STRATEGY.RETAIN_MOST_SUPPORTED;
		Assert.assertEquals(clp.doWork(), 0);
		TestUtils.assertSamFilesSame(clp.OUTPUT, new File(TEST_DATA_DIR, "MarkChimericReads.umiReuse10x.sam"), false);
		Assert.assertTrue(TestUtils.testFilesSame(clp.OUTPUT_REPORT, new File(TEST_DATA_DIR, "umiReuse10x.chimeric_report.txt")));
		Assert.assertTrue(TestUtils.testFilesSame(clp.METRICS, new File(TEST_DATA_DIR, "umiReuse10x.chimeric_read_metrics")));
	}

	@Test(enabled = true)
	public void testPolyT() {
		final MarkChimericReads clp = new MarkChimericReads();
		clp.OUTPUT = TestUtils.getTempReportFile("MarkChimericReads.", ".sam");
		clp.OUTPUT.deleteOnExit();
		clp.OUTPUT_REPORT = TestUtils.getTempReportFile("chimeric_report.", ".txt");
		clp.OUTPUT_REPORT.deleteOnExit();
		clp.METRICS = TestUtils.getTempReportFile("MarkChimericReads.", ".chimeric_read_metrics");
		clp.METRICS.deleteOnExit();
		clp.INPUT = Arrays.asList(new PicardHtsPath(new File(TEST_DATA_DIR, "MarkChimericReads.input.sam")));
		clp.CELL_BC_FILE = new File(TEST_DATA_DIR, "selectedCellBarcodes.txt");
		clp.MARK_UMI_REUSE = false;
		clp.T_RICH_THRESHOLD = 6;
		Assert.assertEquals(clp.doWork(), 0);
		TestUtils.assertSamFilesSame(clp.OUTPUT, new File(TEST_DATA_DIR, "MarkChimericReads.6T.sam"), false);
		Assert.assertTrue(TestUtils.testFilesSame(clp.OUTPUT_REPORT, new File(TEST_DATA_DIR, "6T.chimeric_report.txt")));
		Assert.assertTrue(TestUtils.testFilesSame(clp.METRICS, new File(TEST_DATA_DIR, "6T.chimeric_read_metrics")));
	}

	@Test(enabled = true)
	public void testBoth() {
		final MarkChimericReads clp = new MarkChimericReads();
		clp.OUTPUT = TestUtils.getTempReportFile("MarkChimericReads.", ".sam");
		clp.OUTPUT.deleteOnExit();
		clp.OUTPUT_REPORT = TestUtils.getTempReportFile("chimeric_report.", ".txt");
		clp.OUTPUT_REPORT.deleteOnExit();
		clp.METRICS = TestUtils.getTempReportFile("MarkChimericReads.", ".chimeric_read_metrics");
		clp.METRICS.deleteOnExit();
		clp.INPUT = Arrays.asList(new PicardHtsPath(new File(TEST_DATA_DIR, "MarkChimericReads.input.sam")));
		clp.CELL_BC_FILE = new File(TEST_DATA_DIR, "selectedCellBarcodes.txt");
		clp.T_RICH_THRESHOLD = 6;
		Assert.assertEquals(clp.doWork(), 0);
		TestUtils.assertSamFilesSame(clp.OUTPUT, new File(TEST_DATA_DIR, "MarkChimericReads.both.sam"), false);
		Assert.assertTrue(TestUtils.testFilesSame(clp.OUTPUT_REPORT, new File(TEST_DATA_DIR, "both.chimeric_report.txt")));
		Assert.assertTrue(TestUtils.testFilesSame(clp.METRICS, new File(TEST_DATA_DIR, "both.chimeric_read_metrics")));
	}

}
