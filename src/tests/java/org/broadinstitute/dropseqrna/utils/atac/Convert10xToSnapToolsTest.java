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

package org.broadinstitute.dropseqrna.utils.atac;

import java.io.File;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

public class Convert10xToSnapToolsTest {

  private static final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/atac");
  private static final File TENX_R1_FASTQ = new File(TEST_DATA_DIR, "10X_R1.fastq.gz");
  private static final File TENX_R2_FASTQ = new File(TEST_DATA_DIR, "10X_R2.fastq.gz");
  private static final File SNAPTOOLS_FASTQ = new File(TEST_DATA_DIR, "SnapTools.fastq.gz");
  private static final File ATAC_BARCODES = new File(TEST_DATA_DIR, "10X_atac_barcodes.txt.gz");
  private static final File GEX_BARCODES = new File(TEST_DATA_DIR, "10X_gex_barcodes.txt.gz");
  
  @Test
  public void testDoWork() {
    final Convert10xToSnapTools convert10xToSnapTools = new Convert10xToSnapTools();
    convert10xToSnapTools.INPUT = TENX_R1_FASTQ;
    convert10xToSnapTools.BARCODE_FASTQ = TENX_R2_FASTQ;
    convert10xToSnapTools.BARCODE_ATAC = ATAC_BARCODES;
    convert10xToSnapTools.BARCODE_GEX = GEX_BARCODES;
    convert10xToSnapTools.OUTPUT =
        TestUtils.getTempReportFile("Convert10xToSnapToolsTest.", ".fastq.gz");
    Assert.assertEquals(convert10xToSnapTools.doWork(), 0);
    Assert.assertTrue(TestUtils.testFilesSame(SNAPTOOLS_FASTQ, convert10xToSnapTools.OUTPUT));
  }
}
