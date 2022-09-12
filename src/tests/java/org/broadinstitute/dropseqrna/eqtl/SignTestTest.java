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
package org.broadinstitute.dropseqrna.eqtl;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloseableIterator;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.util.HashMap;
import java.util.Map;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.util.TabbedTextFileWithHeaderParser;
import picard.util.TabbedTextFileWithHeaderParser.Row;

public class SignTestTest {

  private static final File TEST_DATA_DIR =
      new File("testdata/org/broadinstitute/dropseq/eqtl/signtest");

  private static final File EQTL_FILE =
      new File(TEST_DATA_DIR, "input_eqtls.tsv");

  private static final File UNFILTERED_EQTL_FILE =
      new File(TEST_DATA_DIR, "unfiltered_eqtls.tsv");

  @DataProvider
  public Object[][] signTestDataProvider() {
    return new Object[][]{
        {true, true, 0.05, "flip_strip_0.05"},
        {true, false, 0.05, "flip_0.05"},
        {false, true, 0.05, "strip_0.05"},
        {false, false, 0.05, "0.05"},
        {true, true, null, "flip_strip_allq"},
        {true, false, null, "flip_allq"},
        {false, true, null, "strip_allq"},
        {false, false, null, "allq"},
        {false, false, 0.002, "lowq"},
    };
  }

  @Test(dataProvider = "signTestDataProvider")
  public void testSignTest(
      final boolean flipAlleles,
      final boolean stripEnsgSuffix,
      final Double qvalueThreshold,
      final String expectedFilePart
  ) throws IOException {
    final File actualOutputFile =
        File.createTempFile("sign_test_" + expectedFilePart + ".", ".tsv");
    final File expectedOutputFile =
        new File(TEST_DATA_DIR, "sign_test_" + expectedFilePart + ".tsv");

    final SignTest signTest = new SignTest();
    signTest.INPUT = EQTL_FILE;
    signTest.UNFILTERED_EQTL_FILE = UNFILTERED_EQTL_FILE;
    signTest.FLIP_ALLELES = flipAlleles;
    signTest.STRIP_ENSG_SUFFIX = stripEnsgSuffix;
    signTest.QVALUE_THRESHOLD = qvalueThreshold;
    signTest.OUTPUT = actualOutputFile;
    signTest.OUTPUT.deleteOnExit();
    Assert.assertEquals(signTest.doWork(), 0);
    assertMetricsMatch(actualOutputFile, expectedOutputFile);
  }

  private void assertMetricsMatch(final File actualOutputFile, final File expectedOutputFile) {
    final TabbedTextFileWithHeaderParser expectedParser =
        new TabbedTextFileWithHeaderParser(actualOutputFile);
    final TabbedTextFileWithHeaderParser actualParser =
        new TabbedTextFileWithHeaderParser(expectedOutputFile);
    Assert.assertEquals(actualParser.columnLabels(), expectedParser.columnLabels());
    final CloseableIterator<Row> expectedIterator = expectedParser.iterator();
    final CloseableIterator<Row> actualIterator = actualParser.iterator();
    final Row expectedRow = expectedIterator.next();
    final Row actualRow = actualIterator.next();
    Assert.assertFalse(actualIterator.hasNext());
    final Map<String, String> expectedValues = new HashMap<>();
    final Map<String, String> actualValues = new HashMap<>();
    for (final String columnLabel : expectedParser.columnLabels()) {
      if (columnLabel.startsWith("FILE_")) {
        expectedValues.put(columnLabel, new File(expectedRow.getField(columnLabel)).getName());
        actualValues.put(columnLabel, new File(actualRow.getField(columnLabel)).getName());
      } else {
        expectedValues.put(columnLabel, expectedRow.getField(columnLabel));
        actualValues.put(columnLabel, actualRow.getField(columnLabel));
      }
    }
    expectedParser.close();
    actualParser.close();
    Assert.assertEquals(actualValues, expectedValues);
  }
}
