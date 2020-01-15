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
import org.broadinstitute.dropseqrna.barnyard.DigitalExpressionTest;
import org.broadinstitute.dropseqrna.barnyard.SelectCellsByNumTranscripts;
import org.broadinstitute.dropseqrna.barnyard.SingleCellRnaSeqMetricsCollector;
import org.broadinstitute.dropseqrna.utils.SplitBamByCellTest;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.junit.Assert;
import org.testng.annotations.Test;


public class MergeSingleCellRnaSeqMetricsTest {
    private static final File TEST_BAM = SplitBamByCellTest.TEST_BAM;
    private static final File ANNOTATIONS_FILE = new File("/broad/mccarroll/software/metadata/individual_reference/hg19/hg19.refFlat");
    private static final File RIBOSOMAL_INTERVALS_FILE = new File("/broad/mccarroll/software/metadata/individual_reference/hg19/hg19.rRNA.intervals");

    private static File EXPECTED_METRICS = new File ("testdata/org/broadinstitute/dropseq/metrics/MergeSingleCellRnaSeqMetrics.txt");

    @Test
    public void testDoWork() {
        File cellBarcodeFile = TestUtils.getTempReportFile("MergeSingleCellRnaSeqMetrics", ".txt.gz");
        File metrics1 = TestUtils.getTempReportFile("MergeSingleCellRnaSeqMetrics", ".metrics");
        File metrics2 = TestUtils.getTempReportFile("MergeSingleCellRnaSeqMetrics", ".metrics");
        File mergedMetrics = TestUtils.getTempReportFile("MergeSingleCellRnaSeqMetrics", ".metrics");

        List<File> splitBAMFileList = TestUtils.splitBamFile(TEST_BAM, 2);
        Assert.assertEquals(splitBAMFileList.size(), 2);

        selectCellBarcodes(TEST_BAM, cellBarcodeFile);

        // Collect SingleCellRnaSeqMetrics for each of the split BAM files
        collectMetrics(splitBAMFileList.get(0), cellBarcodeFile, metrics1);
        collectMetrics(splitBAMFileList.get(1), cellBarcodeFile, metrics2);

        // Merge SingleCellRnaSeqMetrics
        MergeSingleCellRnaSeqMetrics metricsMerger = new MergeSingleCellRnaSeqMetrics();
        metricsMerger.INPUT = Arrays.asList(metrics1, metrics2);
        metricsMerger.OUTPUT = mergedMetrics;
        Assert.assertEquals(metricsMerger.doWork(), 0);

        Assert.assertTrue(MetricsFile.areMetricsEqual(mergedMetrics, EXPECTED_METRICS));
    }

    private void selectCellBarcodes(final File inputBAM, final File cellListFile) {
        SelectCellsByNumTranscripts cellSelector = new SelectCellsByNumTranscripts();
        String[] args = new String[] {
                "INPUT=" + inputBAM.getAbsolutePath(),
                "MIN_TRANSCRIPTS_PER_CELL=1",
                "OUTPUT=" + cellListFile.getAbsolutePath()
        };

        Assert.assertEquals(cellSelector.instanceMain(args), 0);
    }

    private void collectMetrics(final File inputBAM, final File cellBarcodeFile, final File metricsFile) {
        SingleCellRnaSeqMetricsCollector metricsCollector = new SingleCellRnaSeqMetricsCollector();
        String[] args = new String[] {
                "INPUT=" + inputBAM.getAbsolutePath(),
                "CELL_BARCODE_TAG=XC",
                "CELL_BC_FILE=" + cellBarcodeFile,
                "ANNOTATIONS_FILE=" + ANNOTATIONS_FILE.getAbsolutePath(),
                "RIBOSOMAL_INTERVALS=" + RIBOSOMAL_INTERVALS_FILE.getAbsolutePath(),
                "OUTPUT=" + metricsFile.getAbsolutePath()
        };

        Assert.assertEquals(metricsCollector.instanceMain(args), 0);
    }
}
