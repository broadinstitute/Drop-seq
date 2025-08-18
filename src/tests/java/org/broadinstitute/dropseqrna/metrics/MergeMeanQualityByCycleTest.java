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

import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.SkipException;
import org.testng.annotations.Test;
import picard.analysis.MeanQualityByCycle;
import picard.cmdline.CommandLineProgram;

import java.io.File;
import java.util.Arrays;
import java.util.List;


public class MergeMeanQualityByCycleTest {
    private static final File TEST_BAM = new File("testdata/org/broadinstitute/dropseq/barnyard/DgeStrandFuncTest/DgeStrandFuncTest.bam");
    private static final File EXPECTED_METRICS = new File("testdata/org/broadinstitute/dropseq/metrics/MergeMeanQualityByCycle.txt");

    @Test
    public void testDoWork() {
        if (!CommandLineProgram.checkRInstallation(true)) {
            throw new SkipException(
                    "R is not installed on this machine. It is required for creating the required MeanQualityByCycle chart and running the test."
            );
        }
        File metrics1 = TestUtils.getTempReportFile("MergeMeanQualityByCycle", ".metrics");
        File metrics2 = TestUtils.getTempReportFile("MergeMeanQualityByCycle", ".metrics");
        File mergedMetrics = TestUtils.getTempReportFile("MergeMeanQualityByCycle", ".metrics");

        List<File> splitBAMFileList = TestUtils.splitBamFile(TEST_BAM, 2);
        Assert.assertEquals(splitBAMFileList.size(), 2);

        // Collect MergeMeanQualityByCycle for each of the split BAM files
        collectMetrics(splitBAMFileList.get(0), metrics1);
        collectMetrics(splitBAMFileList.get(1), metrics2);

        // MergeMeanQualityByCycle
        MergeMeanQualityByCycle metricsMerger = new MergeMeanQualityByCycle();
        metricsMerger.INPUT = Arrays.asList(metrics1, metrics2);
        metricsMerger.OUTPUT = mergedMetrics;
        metricsMerger.CHART_OUTPUT = TestUtils.getTempReportFile("MergeMeanQualityByCycle", ".pdf");
        Assert.assertEquals(metricsMerger.doWork(), 0);

        Assert.assertTrue(MetricsFile.areMetricsAndHistogramsEqual(mergedMetrics, EXPECTED_METRICS));
    }

    private void collectMetrics(File inputBAM, File metricsFile) {
        File chartFile = TestUtils.getTempReportFile("MergeMeanQualityByCycle", ".pdf");

        MeanQualityByCycle metricsCollector = new MeanQualityByCycle();
        String[] args = new String[] {
                "INPUT=" + inputBAM.getAbsolutePath(),
                "OUTPUT=" + metricsFile.getAbsolutePath(),
                "CHART_OUTPUT=" + chartFile.getAbsolutePath()
        };
        Assert.assertEquals(metricsCollector.instanceMain(args), 0);
    }
}
