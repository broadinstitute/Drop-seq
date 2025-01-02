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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.dropseqrna.utils.BaseRange;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.nio.PicardHtsPath;

import java.io.File;
import java.util.Collections;
import java.util.List;

public class CorrectAndSplitScrnaReadPairsTest {

 private static final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/beadsynthesis/CorrectAndSplitScrnaReadPairsTest");
 private static final File INPUT_SAM = new File(TEST_DATA_DIR, "correct_barcodes_test.sam");
 private static final File EXPECTED_BARCODES_HIST = new File(TEST_DATA_DIR, "correct_barcodes_test.expected_barcode_metrics.gz");
 private static final String BARCODE_QUALS_TAG = "CY";
 private static final String RAW_BARCODE_TAG = "CR";
 private static final String BASE_RANGE = "1-16";

 @Test(dataProvider = "testBasicDataProvider")
 public void testBasic(boolean tagBothReads) {
  CorrectAndSplitScrnaReadPairs clp = initClp(tagBothReads);
  final File outputSam = new File(clp.OUTPUT_LIST.getParentFile(), "test.0.sam");
  outputSam.deleteOnExit();
  Assert.assertEquals(clp.doWork(), 0);
  final SamReader reader = SamReaderFactory.makeDefault().open(outputSam);
  for (final SAMRecord rec : reader) {
   if (rec.getSecondOfPairFlag() || tagBothReads) {
    final String cellBarcode = rec.getStringAttribute(clp.BARCODE_TAG);
    final String expectedCellBarcode = rec.getReadName().split(":")[1];
    Assert.assertEquals(cellBarcode, expectedCellBarcode, rec.getSAMString());
   }
  }
 }

 @DataProvider(name = "testBasicDataProvider")
 public Object[][] testBasicDataProvider() {
  return new Object[][]{
          {false},
          {true}
  };
 }

 private CorrectAndSplitScrnaReadPairs initClp(boolean tagBothReads) {
  final CorrectAndSplitScrnaReadPairs clp = new CorrectAndSplitScrnaReadPairs();
  clp.INPUT = Collections.singletonList(new PicardHtsPath(INPUT_SAM));
  clp.ALLOWED_BARCODE_COUNTS = EXPECTED_BARCODES_HIST;
  clp.BARCODED_READ = 1;
  clp.BASE_RANGE = BASE_RANGE;
  clp.METRICS = TestUtils.getTempReportFile("test.", ".corrected_barcode_metrics");
  clp.NUM_OUTPUTS = 1;
  clp.OUTPUT_LIST = TestUtils.getTempReportFile("CorrectAndSplitScrnaReadPairsTest.", ".bam_list");
  clp.OUTPUT = new File("test." + clp.OUTPUT_SLUG + ".sam");
  clp.OVERWRITE_EXISTING = true;
  clp.OUTPUT_MANIFEST = TestUtils.getTempReportFile("CorrectAndSplitScrnaReadPairsTest.", ".split_bam_manifest.gz");
  clp.REPORT = TestUtils.getTempReportFile("CorrectAndSplitScrnaReadPairsTest.", ".split_bam_report");
  clp.TAG_BOTH_READS = tagBothReads;
  return clp;
 }

 @Test(dataProvider = "testBasicDataProvider")
 public void testOptionalTags(boolean tagBothReads) {
  CorrectAndSplitScrnaReadPairs clp = initClp(tagBothReads);
  final File outputSam = new File(clp.OUTPUT_LIST.getParentFile(), "test.0.sam");
  outputSam.deleteOnExit();
  clp.BARCODE_QUALS_TAG = BARCODE_QUALS_TAG;
  clp.RAW_BARCODE_TAG = RAW_BARCODE_TAG;
  Assert.assertEquals(clp.doWork(), 0);
  final SamReader reader = SamReaderFactory.makeDefault().open(outputSam);
  final List<BaseRange> baseRanges = org.broadinstitute.dropseqrna.utils.BaseRange.parseBaseRange(BASE_RANGE);

  String rawBarcode = null;
  String barcodeQuals = null;
  for (final SAMRecord rec : reader) {
   if (rec.getFirstOfPairFlag()) {
    rawBarcode = BaseRange.getSequenceForBaseRange(baseRanges, rec.getReadString());
    barcodeQuals = BaseRange.getSequenceForBaseRange(baseRanges, rec.getBaseQualityString());
   }
   if (rec.getSecondOfPairFlag() || tagBothReads) {
    Assert.assertEquals(rec.getAttribute(RAW_BARCODE_TAG), rawBarcode);
    Assert.assertEquals(rec.getAttribute(BARCODE_QUALS_TAG), barcodeQuals);
   }
  }
 }
}
