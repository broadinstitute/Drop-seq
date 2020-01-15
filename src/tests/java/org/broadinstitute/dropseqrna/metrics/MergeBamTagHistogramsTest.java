/*
 * MIT License
 *
 * Copyright 2020 Broad Institute
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

import java.io.File;
import java.util.Arrays;
import java.util.List;

import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.SplitBamByCellTest;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.junit.Assert;
import org.testng.annotations.Test;


public class MergeBamTagHistogramsTest {
    private static final File TEST_BAM = SplitBamByCellTest.TEST_BAM;

    @Test
    public void testDoWork() {
        File originalMetrics = TestUtils.getTempReportFile("MergeBamTagHistograms", ".metrics");
        File metrics1 = TestUtils.getTempReportFile("MergeBamTagHistograms", ".metrics");
        File metrics2 = TestUtils.getTempReportFile("MergeBamTagHistograms", ".metrics");
        File mergedMetrics = TestUtils.getTempReportFile("MergeBamTagHistograms", ".metrics");

        List<File> splitBAMFileList = TestUtils.splitBamFile(TEST_BAM, 2);
        Assert.assertEquals(splitBAMFileList.size(), 2);

        // Compute BamTagHistograms for TEST_BAM and each of the split BAM files
        generateHistograms(TEST_BAM, originalMetrics);
        generateHistograms(splitBAMFileList.get(0), metrics1);
        generateHistograms(splitBAMFileList.get(1), metrics2);

        // Merge BamTagHistograms
        MergeBamTagHistograms histogramMerger = new MergeBamTagHistograms();
        histogramMerger.INPUT = Arrays.asList(metrics1, metrics2);
        histogramMerger.OUTPUT = mergedMetrics;
        Assert.assertEquals(histogramMerger.doWork(), 0);

        ObjectCounter<String> originalCounter = ObjectCounter.readReportFile(originalMetrics);
        ObjectCounter<String> mergedCounter = ObjectCounter.readReportFile(mergedMetrics);
        Assert.assertTrue(mergedCounter.countMapsEqual(originalCounter));
    }

    private void generateHistograms(final File inputBAM, final File metricsFile) {
        BamTagHistogram histogramGenerator = new BamTagHistogram();
        histogramGenerator.INPUT = inputBAM;
        histogramGenerator.TAG = "ZP";
        histogramGenerator.OUTPUT = metricsFile;

        Assert.assertEquals(histogramGenerator.doWork(), 0);
    }
}
