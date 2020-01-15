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

import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.dropseqrna.utils.SplitBamByCellTest;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.junit.Assert;
import org.testng.annotations.Test;
import picard.analysis.CollectRnaSeqMetrics;


public class MergeRnaSeqMetricsTest {
    private static final File TEST_BAM = SplitBamByCellTest.TEST_BAM;
    private static final File REF_FLAT_FILE = new File("/broad/mccarroll/software/metadata/individual_reference/hg19/hg19.refFlat");
    private static final File RIBOSOMAL_INTERVALS_FILE = new File("/broad/mccarroll/software/metadata/individual_reference/hg19/hg19.rRNA.intervals");

    private static File EXPECTED_METRICS = new File ("testdata/org/broadinstitute/dropseq/metrics/MergeRnaSeqMetrics.txt");

    @Test
    public void testDoWork() {
        File metrics1 = TestUtils.getTempReportFile("MergeRnaSeqMetrics", ".metrics");
        File metrics2 = TestUtils.getTempReportFile("MergeRnaSeqMetrics", ".metrics");
        File mergedMetrics = TestUtils.getTempReportFile("MergeRnaSeqMetrics", ".metrics");

        List<File> splitBAMFileList = TestUtils.splitBamFile(TEST_BAM, 2);
        Assert.assertEquals(splitBAMFileList.size(), 2);

        // Run CollectRnaSeqMetrics on each of the split BAM files
        collectMetrics(splitBAMFileList.get(0), metrics1);
        collectMetrics(splitBAMFileList.get(1), metrics2);

        // Merge RnaSeq mettrics
        MergeRnaSeqMetrics metricsMerger = new MergeRnaSeqMetrics();
        metricsMerger.INPUT = Arrays.asList(metrics1, metrics2);
        metricsMerger.OUTPUT = mergedMetrics;
        Assert.assertEquals(metricsMerger.doWork(), 0);

        Assert.assertTrue(MetricsFile.areMetricsEqual(mergedMetrics, EXPECTED_METRICS));
    }

    private void collectMetrics(File inputBAM, File metricsFile) {
        CollectRnaSeqMetrics metricsCollector = new CollectRnaSeqMetrics();
        String[] args = new String[] {
                "INPUT=" + inputBAM.getAbsolutePath(),
                "REF_FLAT=" + REF_FLAT_FILE.getAbsolutePath(),
                "RIBOSOMAL_INTERVALS=" + RIBOSOMAL_INTERVALS_FILE.getAbsolutePath(),
                "STRAND_SPECIFICITY=NONE",
                "OUTPUT=" + metricsFile.getAbsolutePath()
        };
        Assert.assertEquals(metricsCollector.instanceMain(args), 0);
    }
}
