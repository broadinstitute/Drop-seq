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
package org.broadinstitute.dropseqrna.beadsynthesis;

import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.Collections;

public class CountBarcodeSequencesTest {
 private static final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/beadsynthesis/CountBarcodeSequencesTest");
 private static final File PAIRED_BAM = new File(TEST_DATA_DIR, "unmapped_paired_reads.sam");
 private static final File ALLOW_LIST = new File(TEST_DATA_DIR, "allowList.txt");
 private static final File EXPECTED_NO_ALLOW_LIST = new File(TEST_DATA_DIR, "noAllowList.count_barcode_sequences_metrics");
 private static final File EXPECTED_ALLOW_LIST = new File(TEST_DATA_DIR, "allowList.count_barcode_sequences_metrics");
 private static final File EXPECTED_ALLOW_LIST_PSEUDOCOUNT = new File(TEST_DATA_DIR, "allowListPseudoCount.count_barcode_sequences_metrics");

 @Test(dataProvider = "testBasicDataProvider")
 public void testBasic(final String testName, final File allowList, final int allowListPseudoCount, final File expectedMetricsFile, final int numThreads)
         throws FileNotFoundException {
  final CountBarcodeSequences clp = new CountBarcodeSequences();
  clp.BARCODED_READ = 1;
  clp.OUTPUT = TestUtils.getTempReportFile(testName + ".", ".count_barcode_sequences_metrics");
  clp.INPUT = Collections.singletonList(PAIRED_BAM);
  clp.NUM_THREADS = numThreads;
  clp.BASE_RANGE = "1-4";
  clp.ALLOWED_BARCODES = allowList;
  clp.ALLOWLIST_PSEUDOCOUNT = allowListPseudoCount;
  Assert.assertEquals(clp.doWork(), 0);
  Assert.assertTrue(TestUtils.testMetricsFilesEqual(expectedMetricsFile, clp.OUTPUT));
 }

 @DataProvider(name = "testBasicDataProvider")
 public Object[][] testBasicDataProvider() {
  return new Object[][] {
          {"noAllowList", null, 0, EXPECTED_NO_ALLOW_LIST, 1},
          {"allowList", ALLOW_LIST, 0, EXPECTED_ALLOW_LIST, 1},
          {"allowListPseudoCount", ALLOW_LIST, 1, EXPECTED_ALLOW_LIST_PSEUDOCOUNT, 1},
          {"noAllowListThreads", null, 0, EXPECTED_NO_ALLOW_LIST, 4},
          {"allowListThreads", ALLOW_LIST, 0, EXPECTED_ALLOW_LIST, 4},
          {"allowListPseudoCountThreads", ALLOW_LIST, 1, EXPECTED_ALLOW_LIST_PSEUDOCOUNT, 4},
  };

 }

 @Test(dataProvider = "twoInputsDataProvider")
 public void testTwoInputsThreadCountsAgree(final File allowList) throws FileNotFoundException {
  // {NUM_THREADS, READER_THREADS, COUNTING_THREADS}: the first is the single-threaded reference.  NUM_THREADS is
  // ignored when the other two are given, and the split is automatic when they are null.
  final Integer[][] threads = {{1, null, null}, {2, null, null}, {3, null, null}, {4, null, null}, {8, null, null},
          {1, 2, 0}, {1, 1, 1}, {1, 3, 1}, {1, 2, 4}};
  File reference = null;
  for (final Integer[] t : threads) {
   final CountBarcodeSequences clp = new CountBarcodeSequences();
   clp.BARCODED_READ = 1;
   clp.OUTPUT = TestUtils.getTempReportFile("twoInputs" + t[0] + "." + t[1] + "." + t[2] + ".", ".count_barcode_sequences_metrics");
   clp.INPUT = Arrays.asList(PAIRED_BAM, PAIRED_BAM);
   clp.NUM_THREADS = t[0];
   clp.READER_THREADS = t[1];
   clp.COUNTING_THREADS = t[2];
   clp.ALLOWED_BARCODES = allowList;
   clp.BASE_RANGE = "1-4";
   Assert.assertEquals(clp.doWork(), 0);
   if (reference == null) {
    reference = clp.OUTPUT;
   } else {
    Assert.assertTrue(TestUtils.testMetricsFilesEqual(reference, clp.OUTPUT));
   }
  }
 }

 @DataProvider(name = "twoInputsDataProvider")
 public Object[][] twoInputsDataProvider() {
  return new Object[][] {{null}, {ALLOW_LIST}};
 }

 @Test
 public void testInvalidNumThreads() {
  final CountBarcodeSequences clp = new CountBarcodeSequences();
  clp.BARCODED_READ = 1;
  clp.OUTPUT = TestUtils.getTempReportFile("badThreads.", ".count_barcode_sequences_metrics");
  clp.INPUT = Collections.singletonList(PAIRED_BAM);
  clp.BASE_RANGE = "1-4";
  clp.NUM_THREADS = 0;
  final String[] errors = clp.customCommandLineValidation();
  Assert.assertNotNull(errors);
  Assert.assertTrue(Arrays.stream(errors).anyMatch(e -> e.contains("NUM_THREADS")));
 }

 @Test(dataProvider = "invalidThreadsDataProvider")
 public void testInvalidReaderAndCountingThreads(final Integer readerThreads, final Integer countingThreads) {
  final CountBarcodeSequences clp = new CountBarcodeSequences();
  clp.BARCODED_READ = 1;
  clp.OUTPUT = TestUtils.getTempReportFile("badReaderCountingThreads.", ".count_barcode_sequences_metrics");
  clp.INPUT = Collections.singletonList(PAIRED_BAM);
  clp.BASE_RANGE = "1-4";
  clp.READER_THREADS = readerThreads;
  clp.COUNTING_THREADS = countingThreads;
  final String[] errors = clp.customCommandLineValidation();
  Assert.assertNotNull(errors);
  Assert.assertTrue(Arrays.stream(errors).anyMatch(e -> e.contains("READER_THREADS")));
 }

 @DataProvider(name = "invalidThreadsDataProvider")
 public Object[][] invalidThreadsDataProvider() {
  return new Object[][] {{2, null}, {null, 1}, {0, 1}, {2, -1}};
 }
}
