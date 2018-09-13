/*
 * MIT License
 *
 * Copyright 2018 Broad Institute
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
package org.broadinstitute.dropseqrna.utils.editdistance;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import picard.util.TabbedTextFileWithHeaderParser;

public class DetectBeadSubstitutionErrorsTest {
    private static final File TEST_FILE = new File("testdata/org/broadinstitute/transcriptome/utils/editdistance/DetectBeadSubstitutionErrors.bam");

    private static final String BIG_BARCODE = "AGTGAGACAAGG";
    private static final Set<String> SMALL_BARCODES = new HashSet<>(Arrays.asList(
            "AGTGCGACAAGG",
            "ACTGAGACAAGG",
            "GGTGAGACAAGG"
    ));

    /**
     * Very lame test -- just confirms that program doesn't crash and return 0 exit status.
     */
    @Test
    public void testBasic() throws IOException {
        final DetectBeadSubstitutionErrors clp = new DetectBeadSubstitutionErrors();
        clp.INPUT = Arrays.asList(TEST_FILE);
        clp.OUTPUT_REPORT = File.createTempFile("DetectBeadSubstitutionErrorsTest.", ".substitution_report.txt");
        clp.OUTPUT = File.createTempFile("DetectBeadSubstitutionErrorsTest.", ".sam");
        clp.OUTPUT_REPORT.deleteOnExit();
        clp.OUTPUT_SUMMARY = File.createTempFile("DetectBeadSubstitutionErrorsTest.", "substitution_summary.txt");
        clp.OUTPUT.deleteOnExit();
        Assert.assertEquals(clp.doWork(), 0);
        final TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(clp.OUTPUT_REPORT);
        int numRows = 0;
        for (final TabbedTextFileWithHeaderParser.Row row : parser) {
            ++numRows;
            Assert.assertEquals(row.getField("intended_barcode"), BIG_BARCODE);
            final String neighbor_barcode = row.getField("neighbor_barcode");
            Assert.assertTrue(SMALL_BARCODES.contains(neighbor_barcode));
            Assert.assertTrue(Boolean.parseBoolean(row.getField("repaired")));
        }
        final SamReader reader = SamReaderFactory.makeDefault().open(clp.OUTPUT);
        for (final SAMRecord rec : reader)
			Assert.assertEquals(rec.getStringAttribute("XC"), BIG_BARCODE);
        CloserUtil.close(reader);
    }
}
