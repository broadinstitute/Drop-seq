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
package org.broadinstitute.dropseqrna.utils;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;


public class MergeBaseDistributionAtReadPositionTest {
    private static final File TEST_BAM = new File("testdata/org/broadinstitute/dropseq/barnyard/DgeStrandFuncTest/DgeStrandFuncTest.bam");

    @Test
    public void testDoWork() {
        File originalMetrics = TestUtils.getTempReportFile("MergeBaseDistributionAtReadPosition", ".metrics");
        File metrics1 = TestUtils.getTempReportFile("MergeBaseDistributionAtReadPosition", ".metrics");
        File metrics2 = TestUtils.getTempReportFile("MergeBaseDistributionAtReadPosition", ".metrics");
        File mergedMetrics = TestUtils.getTempReportFile("MergeBaseDistributionAtReadPosition", ".metrics");

        List<File> splitBAMFileList = TestUtils.splitBamFile(TEST_BAM, 2);
        Assert.assertEquals(splitBAMFileList.size(), 2);

        // Compute BaseDistributionAtReadPosition for TEST_BAM and each of the split BAM files
        computeDistribution(TEST_BAM, originalMetrics);
        computeDistribution(splitBAMFileList.get(0), metrics1);
        computeDistribution(splitBAMFileList.get(1), metrics2);

        // Merge BaseDistributionAtReadPosition
        MergeBaseDistributionAtReadPosition merger = new MergeBaseDistributionAtReadPosition();
        merger.INPUT = Arrays.asList(metrics1, metrics2);
        merger.OUTPUT = mergedMetrics;
        Assert.assertEquals(merger.doWork(), 0);

        BaseDistributionMetricCollection originalMetricCollection = BaseDistributionMetricCollection.readBaseDistribution(originalMetrics);
        BaseDistributionMetricCollection mergedMetricCollection = BaseDistributionMetricCollection.readBaseDistribution(mergedMetrics);
        boolean sameDistributions = originalMetricCollection.distributionsEqual(mergedMetricCollection);
        Assert.assertTrue(sameDistributions);
    }

    private void computeDistribution(File inputBAM, File metricsFile) {
        BaseDistributionAtReadPosition program = new BaseDistributionAtReadPosition();
        final String args[] = {
                "INPUT=" + inputBAM.getAbsolutePath(),
                "TAG=" + "XC",
                "OUTPUT=" + metricsFile.getAbsolutePath()
        };

        Assert.assertEquals(program.instanceMain(args), 0);
    }
}
