/*
 * MIT License
 *
 * Copyright 2025 Broad Institute
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
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.nio.PicardHtsPath;

import java.util.Collections;

public class CorrectScrnaReadPairsTest {
    @Test(dataProvider = "testBasicDataProvider")
    public void testBasic(boolean tagBothReads) {
        CorrectScrnaReadPairs clp = new CorrectScrnaReadPairs();
        clp.INPUT = Collections.singletonList(new PicardHtsPath(CorrectAndSplitScrnaReadPairsTest.INPUT_SAM));
        clp.ALLOWED_BARCODE_COUNTS = CorrectAndSplitScrnaReadPairsTest.EXPECTED_BARCODES_HIST;
        clp.BARCODED_READ = 1;
        clp.BASE_RANGE = CorrectAndSplitScrnaReadPairsTest.BASE_RANGE;
        clp.METRICS = TestUtils.getTempReportFile("test.", ".corrected_barcode_metrics");
        clp.OUTPUT = TestUtils.getTempReportFile("CorrectScrnaReadPairsTest.", ".sam");
        clp.TAG_BOTH_READS = tagBothReads;
        Assert.assertEquals(clp.doWork(), 0);
        final SamReader reader = SamReaderFactory.makeDefault().open(clp.OUTPUT);
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
}
