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
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Collections;

public class CorrectAndSplitScrnaReadPairsTest {

 private static final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/beadsynthesis/CorrectAndSplitScrnaReadPairsTest");
 private static final File INPUT_SAM = new File(TEST_DATA_DIR, "correct_barcodes_test.sam");
 private static final File EXPECTED_BARCODES_HIST = new File(TEST_DATA_DIR, "correct_barcodes_test.expected_barcode_metrics.gz");

 @Test
 public void testBasic() {
  final CorrectAndSplitScrnaReadPairs clp = new CorrectAndSplitScrnaReadPairs();
  clp.INPUT = Collections.singletonList(INPUT_SAM);
  clp.ALLOWED_BARCODE_COUNTS = EXPECTED_BARCODES_HIST;
  clp.BARCODED_READ = 1;
  clp.BASE_RANGE = "1-16";
  clp.METRICS = TestUtils.getTempReportFile("test.", ".corrected_barcode_metrics");
  clp.NUM_OUTPUTS = 1;
  clp.OUTPUT_LIST = TestUtils.getTempReportFile( "CorrectAndSplitScrnaReadPairsTest.", ".bam_list");
  clp.OUTPUT = new File("test." + clp.OUTPUT_SLUG + ".sam");
  clp.OVERWRITE_EXISTING = true;
  final File outputSam = new File(clp.OUTPUT_LIST.getParentFile(), "test.0.sam");
  outputSam.deleteOnExit();
  clp.OUTPUT_MANIFEST = TestUtils.getTempReportFile("CorrectAndSplitScrnaReadPairsTest.", ".split_bam_manifest.gz");
  clp.REPORT = TestUtils.getTempReportFile("CorrectAndSplitScrnaReadPairsTest.", ".split_bam_report");
  Assert.assertEquals(clp.doWork(), 0);
  final SamReader reader = SamReaderFactory.makeDefault().open(outputSam);
  for (final SAMRecord rec: reader) {
   if (rec.getSecondOfPairFlag()) {
    final String cellBarcode = rec.getStringAttribute(clp.BARCODE_TAG);
    final String expectedCellBarcode = rec.getReadName().split(":")[1];
    Assert.assertEquals(cellBarcode, expectedCellBarcode, rec.getSAMString());
   }
  }
 }
}
