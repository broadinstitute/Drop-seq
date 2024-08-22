/*
 * MIT License
 *
 * Copyright 2022 Broad Institute
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

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IterableAdapter;
import org.broadinstitute.dropseqrna.barnyard.ChimericReportEditDistanceCollapse;
import org.broadinstitute.dropseqrna.barnyard.GatherMolecularBarcodeDistributionByGene;
import org.broadinstitute.dropseqrna.barnyard.MarkChimericReads;
import org.broadinstitute.dropseqrna.barnyard.UMICollectionByCellParser;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.List;

public class ChimericReportEditDistanceCollapseTest {
 private static final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/barnyard");
 private static final File CHIMERIC_REPORT = new File(TEST_DATA_DIR, "umiReuse.chimeric_report.txt");
 private static final File TWO_CELLS_CHIMERIC_REPORT = new File(TEST_DATA_DIR, "two_cells.chimeric_report.txt");
 private static final File ED1_CHIMERIC_REPORT = new File(TEST_DATA_DIR, "ed_collapse.chimeric_report.txt");

 private static final File CELL_BC_FILE = new File(TEST_DATA_DIR, "selectedCellBarcodes.txt");

 @Test(dataProvider = "testReportParserDataProvider")
 public void testReportParser(boolean skipChimeric) {
  final UMICollectionByCellParser parser = new UMICollectionByCellParser(CHIMERIC_REPORT, skipChimeric);
  int numUmis = 0;
  for (final List<UMICollection> umiCollections : new IterableAdapter<>(parser)) {
   for (final UMICollection umiCollection: umiCollections) {
    numUmis += umiCollection.getMolecularBarcodeCounts().getSize();
   }
   CloserUtil.close(parser);
   int expectNumUmis = 0;
   final TabbedTextFileWithHeaderParser reportParser = new TabbedTextFileWithHeaderParser(CHIMERIC_REPORT);
   for (TabbedTextFileWithHeaderParser.Row row : reportParser) {
    if (skipChimeric && row.getField(MarkChimericReads.CHIMERIC_COLUMN).equals("true")) {
     continue;
    }
    ++expectNumUmis;
   }
   CloserUtil.close(reportParser);
   Assert.assertEquals(numUmis, expectNumUmis);
  }
 }

 @DataProvider(name = "testReportParserDataProvider")
 Object[][] testReportParserDataProvider() {
  return new Object[][]{
          {false},
          {true},
  };
 }

 @Test(dataProvider = "testBasicDataProvider")
 public void testBasic(boolean ignoreChimeric, int editDistance) {
  final ChimericReportEditDistanceCollapse clp = new ChimericReportEditDistanceCollapse();
  clp.INPUT = CHIMERIC_REPORT;
  clp.OUTPUT = TestUtils.getTempReportFile("ChimericReportEditDistanceCollapseTest", ".molBC.txt.gz");
  clp.EDIT_DISTANCE = editDistance;
  clp.IGNORE_CHIMERIC = ignoreChimeric;
  Assert.assertEquals(clp.doWork(), 0);
  // There is no edit distance collapse in this input file, so just confirm that number of rows is correct
  // based on ignoreChimeric
  final List<TabbedTextFileWithHeaderParser.Row> inputLines = slurpRows(clp.INPUT);
  final long numChimeric = inputLines.stream().filter(row -> row.getField(MarkChimericReads.CHIMERIC_COLUMN).equals("true")).count();
  final long outputRows = slurpRows(clp.OUTPUT).size();
  final long expectedRows = ignoreChimeric? inputLines.size() - numChimeric: inputLines.size();
  Assert.assertEquals(outputRows, expectedRows);
 }
 @DataProvider(name = "testBasicDataProvider")
 Object[][] testBasicDataProvider() {
  return new Object[][]{
          {false, 0},
          {true, 0},
          {false, 1},
          {true, 1},
  };
 }

 @Test
 public void testSelectedCells() throws FileNotFoundException {
  final ChimericReportEditDistanceCollapse clp = new ChimericReportEditDistanceCollapse();
  clp.INPUT = TWO_CELLS_CHIMERIC_REPORT;
  clp.CELL_BC_FILE = CELL_BC_FILE;
  clp.IGNORE_CHIMERIC = false;
  // not gzipped so it can be read by slurpLines
  clp.OUTPUT = TestUtils.getTempReportFile("ChimericReportEditDistanceCollapseTest", ".molBC.txt");
  Assert.assertEquals(clp.doWork(), 0);
  final List<String> reportLines = IOUtil.slurpLines(clp.OUTPUT);
  // header line plus one of the two cells
  Assert.assertEquals(reportLines.size(), 2);
 }

 @Test(dataProvider = "testEditDistanceCollapse")
 public void testEditDistanceCollapse(boolean ignoreChimeric, int expectedObservations) {
  // TODO
  final ChimericReportEditDistanceCollapse clp = new ChimericReportEditDistanceCollapse();
  clp.INPUT = ED1_CHIMERIC_REPORT;
  clp.EDIT_DISTANCE = 1;
  clp.IGNORE_CHIMERIC = ignoreChimeric;
  clp.OUTPUT = TestUtils.getTempReportFile("ChimericReportEditDistanceCollapseTest", ".molBC.txt.gz");
  Assert.assertEquals(clp.doWork(), 0);
  final TabbedTextFileWithHeaderParser.Row row;
  try (CloseableIterator<TabbedTextFileWithHeaderParser.Row> it = new TabbedTextFileWithHeaderParser(clp.OUTPUT).iterator()) {
   row = it.next();
   Assert.assertFalse(it.hasNext());
  }
  Assert.assertEquals(Integer.parseInt(row.getField(GatherMolecularBarcodeDistributionByGene.NUM_OBS_COLUMN)), expectedObservations);

 }

 @DataProvider(name = "testEditDistanceCollapse")
 Object[][] testEditDistanceCollapse() {
  return new Object[][] {
          {true, 18},
          {false, 19},
  };
 }

 private static List<TabbedTextFileWithHeaderParser.Row> slurpRows(final File f) {
  try (final CloseableIterator<TabbedTextFileWithHeaderParser.Row> inputIterator = new TabbedTextFileWithHeaderParser(f).iterator()) {
   return inputIterator.toList();
  }
 }
}
