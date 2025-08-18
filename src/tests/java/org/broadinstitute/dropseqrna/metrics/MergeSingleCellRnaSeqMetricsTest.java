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

import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;


public class MergeSingleCellRnaSeqMetricsTest {
    private static final File SPLIT_BAM_METRICS1 = new File("testdata/org/broadinstitute/dropseq/metrics/MergeSingleCellRnaSeqSplitBamMetrics1.txt");
    private static final File SPLIT_BAM_METRICS2 = new File("testdata/org/broadinstitute/dropseq/metrics/MergeSingleCellRnaSeqSplitBamMetrics2.txt");
    private static final File EXPECTED_METRICS = new File("testdata/org/broadinstitute/dropseq/metrics/MergeSingleCellRnaSeqMetrics.txt");

    @Test
    public void testDoWork() {
        File mergedMetrics = TestUtils.getTempReportFile("MergeSingleCellRnaSeqMetrics", ".metrics");

        // Merge SingleCellRnaSeqMetrics
        MergeSingleCellRnaSeqMetrics metricsMerger = new MergeSingleCellRnaSeqMetrics();
        metricsMerger.INPUT = Arrays.asList(SPLIT_BAM_METRICS1, SPLIT_BAM_METRICS2);
        metricsMerger.OUTPUT = mergedMetrics;
        Assert.assertEquals(metricsMerger.doWork(), 0);

        Assert.assertTrue(MetricsFile.areMetricsEqual(mergedMetrics, EXPECTED_METRICS));
    }
}
