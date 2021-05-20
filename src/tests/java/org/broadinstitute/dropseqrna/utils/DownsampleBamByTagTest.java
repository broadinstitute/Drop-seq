/*
 * MIT License
 *
 * Copyright 2021 Broad Institute
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
package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.dropseqrna.metrics.BamTagHistogram;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class DownsampleBamByTagTest {

    private static final File TESTDATA_DIR = new File("testdata/org/broadinstitute/dropseq/utils");
    private static final File ALIGNED_UNPAIRED_BAM = new File(TESTDATA_DIR, "N701_small.bam");
    private static final String TAG = "XC";
    private static final File ALIGNED_PAIRED_BAM = new File(TESTDATA_DIR, "d0GRIA3_A.multi_organism.MOUSE.census.paired.bam");

    // TODO: paired test

    @Test(dataProvider = "testUnpairedDataProvider")
    public void testUnpaired(boolean USE_PROBABILISTIC_STRATEGY, int READ_MQ, boolean FILTER_PCR_DUPLICATES, boolean use_NUM_READS) {
        final DownsampleBamByTag clp = new DownsampleBamByTag();
        clp.USE_PROBABILISTIC_STRATEGY = USE_PROBABILISTIC_STRATEGY;
        clp.READ_MQ = READ_MQ;
        clp.FILTER_PCR_DUPLICATES = FILTER_PCR_DUPLICATES;
        clp.PAIRED_READS = false;
        clp.TAG = TAG;
        clp.INPUT = ALIGNED_UNPAIRED_BAM;
        clp.OUTPUT = TestUtils.getTempReportFile("testUnpaired.", ".sam");
        ObjectCounter<String> originalCounts = getTagCounts(clp.INPUT, TAG, READ_MQ, FILTER_PCR_DUPLICATES);
        final ObjectCounter<String> expectedCounts = new ObjectCounter<>();
        if (use_NUM_READS) {
            clp.NUM_READS = 10; //originalCounts.getCountForKey(originalCounts.getMin());
            for (final String tagValue: originalCounts.getKeys()) {
                expectedCounts.setCount(tagValue, clp.NUM_READS);
            }
        } else {
            int tagIndex = 5;
            final int numTagsDesired = 10;
            for (final String tagValue: originalCounts.getKeys()) {
                expectedCounts.setCount(tagValue, ++tagIndex);
                if (expectedCounts.getSize() >= numTagsDesired) break;
            }
            clp.TAG_FILE = TestUtils.getTempReportFile("testUnpaired.", ".XC.txt");
            writeTagValuesToFile(clp.TAG_FILE, expectedCounts);
        }
        Assert.assertEquals(clp.doWork(), 0);
        // Quality and dupe filtering alredy done, so shouldn't need to do it here.
        ObjectCounter<String> actualCounts = getTagCounts(clp.OUTPUT, TAG, 0, false);
        final ArrayList<String> expectedKeys = new ArrayList<>(expectedCounts.getKeys());
        Collections.sort(expectedKeys);
        final ArrayList<String> actualKeys = new ArrayList<>(actualCounts.getKeys());
        Collections.sort(actualKeys);
        if (USE_PROBABILISTIC_STRATEGY) {
            // It appears that with probabilistic, and these small numbers, that some barcodes will be lost
            Assert.assertTrue(actualKeys.size() <= expectedKeys.size());

        } else {
            Assert.assertEquals(actualKeys, expectedKeys);
            // Can't check for match if probabilistic
            for (final String tagValue : actualCounts.getKeys()) {
                Assert.assertEquals(actualCounts.getCountForKey(tagValue),
                        Math.min(expectedCounts.getCountForKey(tagValue), originalCounts.getCountForKey(tagValue)));
            }
        }
    }

    @DataProvider(name="testUnpairedDataProvider")
    public Object[][] testUnpairedDataProvider() {
        final List<Object[]> ret = new ArrayList<>();
        final boolean[] tf = {true, false};
        final int[] readMqValues = {10, 0};
        for (final boolean USE_PROBABILISTIC_STRATEGY: tf) {
            for (final int READ_MQ: readMqValues) {
                for (final boolean FILTER_PCR_DUPLICATES : tf) {
                    for (final boolean use_NUM_READS: tf) {
                        final Object[] parameters = {USE_PROBABILISTIC_STRATEGY, READ_MQ, FILTER_PCR_DUPLICATES, use_NUM_READS};
                        ret.add(parameters);
                    }
                }
            }
        }
        return ret.toArray(new Object[0][]);
    }

    public static ObjectCounter<String> getTagCounts(final File bam, final String tag, int READ_MQ, boolean FILTER_PCR_DUPLICATES) {
        return new BamTagHistogram().getBamTagCounts (bam, tag, READ_MQ, FILTER_PCR_DUPLICATES);
    }

    private void writeTagValuesToFile(final File outFile, final ObjectCounter<String> desiredCounts) {
        final PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(outFile));
        for (final String tagValue : desiredCounts.getKeys()) {
            out.println(StringUtil.join("\t", tagValue, desiredCounts.getCountForKey(tagValue)));
        }
        out.close();
    }

    // Very rudimentary test -- just runs the code and confirms it doesn't crash.
    @Test
    public void testPaired() {
        final DownsampleBamByTag clp = new DownsampleBamByTag();
        clp.USE_PROBABILISTIC_STRATEGY = true;
        clp.PAIRED_READS = true;
        clp.TAG = "RG";
        clp.NUM_READS = 10;
        clp.INPUT = ALIGNED_PAIRED_BAM;
        clp.OUTPUT = TestUtils.getTempReportFile("testPaired.", ".sam");
        Assert.assertEquals(clp.doWork(), 0);
    }
}
