/*
 * MIT License
 *
 * Copyright 2023 Broad Institute
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
package org.broadinstitute.dropseqrna.metrics;

import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.nio.PicardHtsPath;

import java.io.File;
import java.util.List;

public class GatherUMIReadIntervalsTest {
 static final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/metrics/GatherUMIReadIntervals");
 // Just a random aligned SAM, and corresponding selected cells
 private static final File INPUT = new File(TEST_DATA_DIR, "gatherUmiReadIntervals.sam");
 private static final File SELECTED_CELLS = new File(TEST_DATA_DIR, "gatherUmiReadIntervals.selectedCellBarcodes.txt");

 static final File EXPECTED_REPORT = new File(TEST_DATA_DIR, "GatherUMIReadIntervalsTest.tsv");
 static final File EXPECTED_REPORT_ED3 = new File(TEST_DATA_DIR, "GatherUMIReadIntervalsTest.ed3.tsv");

 @Test(dataProvider = "testBasicDataProvider")
 public void testBasic(final int editDistance, final File expectedReport) {
  final GatherUMIReadIntervals clp = new GatherUMIReadIntervals();
  clp.INPUT = List.of(new PicardHtsPath(INPUT));
  clp.CELL_BC_FILE = SELECTED_CELLS;
  clp.OUTPUT = TestUtils.getTempReportFile("GatherUMIReadIntervalsTest.", ".tsv");
  clp.EDIT_DISTANCE = editDistance;
  Assert.assertEquals(clp.doWork(), 0);
  Assert.assertTrue(TestUtils.testFilesSame(expectedReport, clp.OUTPUT));
 }

 @DataProvider(name = "testBasicDataProvider")
 Object[][] testBasicDataProvider() {
  return new Object[][] {
          {0, EXPECTED_REPORT},
          {3, EXPECTED_REPORT_ED3},
  };
 }

}
