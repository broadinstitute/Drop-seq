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
package org.broadinstitute.dropseqrna.metrics;

import htsjdk.samtools.util.CloseableIterator;
import org.broadinstitute.dropseqrna.utils.TagBamWithReadSequenceExtended;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.File;
import java.util.Arrays;

public class MergeTagBamWithBarcodeSummariesTest {

    private static final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/utils");
    private static final File SUMMARY = new File(TEST_DATA_DIR, "paired_reads.tagged_Cellular.bam_summary.txt");

    @Test
    public void testBasic() {
        final MergeTagBamWithBarcodeSummaries clp = new MergeTagBamWithBarcodeSummaries();
        clp.INPUT = Arrays.asList(SUMMARY, SUMMARY);
        clp.OUTPUT = TestUtils.getTempReportFile("MergeTagBamWithBarcodeSummariesTest.", ".bam_summary.txt");
        Assert.assertEquals(clp.doWork(), 0);
        try (
                final TabbedTextFileWithHeaderParser mergedParser = new TabbedTextFileWithHeaderParser(clp.OUTPUT);
                final TabbedTextFileWithHeaderParser originalParser = new TabbedTextFileWithHeaderParser(SUMMARY);
                final CloseableIterator<TabbedTextFileWithHeaderParser.Row> mergedIt = mergedParser.iterator();
                final CloseableIterator<TabbedTextFileWithHeaderParser.Row> originalIt = originalParser.iterator()
        ) {
            while (originalIt.hasNext()) {
                Assert.assertTrue(mergedIt.hasNext());
                final TabbedTextFileWithHeaderParser.Row originalRow = originalIt.next();
                final TabbedTextFileWithHeaderParser.Row mergedRow = mergedIt.next();
                Assert.assertEquals(mergedRow.getField(TagBamWithReadSequenceExtended.NUM_FAILED_BASES_COLUMN),
                        originalRow.getField(TagBamWithReadSequenceExtended.NUM_FAILED_BASES_COLUMN));
                Assert.assertEquals(mergedRow.getIntegerField(TagBamWithReadSequenceExtended.NUM_BARCODES_COLUMN).intValue(),
                        2 * originalRow.getIntegerField(TagBamWithReadSequenceExtended.NUM_BARCODES_COLUMN));
            }
        }
    }
}
