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
package org.broadinstitute.dropseqrna.barnyard;

import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.dropseqrna.utils.SplitBamByCellTest;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.junit.Assert;
import org.testng.annotations.Test;
import picard.annotation.LocusFunction;


public class MergeCellsByNumTranscriptsTest {

    private static final File TEST_BAM = SplitBamByCellTest.TEST_BAM;

    @Test
    public void testDoWork() {
        File originalCellBarcodeFile = TestUtils.getTempReportFile("MergeCellsByNumTranscripts", ".txt.gz");
        File cellBarcodeFile1 = TestUtils.getTempReportFile("MergeCellsByNumTranscripts", ".txt.gz");
        File cellBarcodeFile2 = TestUtils.getTempReportFile("MergeCellsByNumTranscripts", ".txt.gz");
        File mergedCellBarcodeFile = TestUtils.getTempReportFile("MergeCellsByNumTranscripts", ".txt.gz");
        File originalMetrics = TestUtils.getTempReportFile("MergeCellsByNumTranscripts", ".metrics");
        File metrics1 = TestUtils.getTempReportFile("MergeCellsByNumTranscripts", ".metrics");
        File metrics2 = TestUtils.getTempReportFile("MergeCellsByNumTranscripts", ".metrics");
        File mergedMetrics = TestUtils.getTempReportFile("MergeCellsByNumTranscripts", ".metrics");

        List<File> splitBAMFileList = TestUtils.splitBamFile(TEST_BAM, 2);
        Assert.assertEquals(splitBAMFileList.size(), 2);

        // Run SelectCellsByNumTranscripts for TEST_BAN as well as each of the split BAM files
        selectCells(TEST_BAM, originalCellBarcodeFile, originalMetrics);
        selectCells(splitBAMFileList.get(0), cellBarcodeFile1, metrics1);
        selectCells(splitBAMFileList.get(1), cellBarcodeFile2, metrics2);

        // MergeCellsByNumTranscripts
        MergeCellsByNumTranscripts cellMerger = new MergeCellsByNumTranscripts();
        cellMerger.INPUT = Arrays.asList(cellBarcodeFile1, cellBarcodeFile2);
        cellMerger.INPUT_METRICS = Arrays.asList(metrics1, metrics2);
        cellMerger.OUTPUT = mergedCellBarcodeFile;
        cellMerger.METRICS = mergedMetrics;
        Assert.assertEquals(cellMerger.doWork(), 0);

        Set<String> originalCellBarcodes = new HashSet<String>(SelectCellsByNumTranscripts.readBarcodes(originalCellBarcodeFile));
        Set<String> mergedCellBarcodes = new HashSet<String>(SelectCellsByNumTranscripts.readBarcodes(mergedCellBarcodeFile));
        Assert.assertTrue(
                mergedCellBarcodes.size() == originalCellBarcodes.size() &&
                originalCellBarcodes.containsAll(mergedCellBarcodes));

        Assert.assertTrue(MetricsFile.areMetricsEqual(mergedMetrics, originalMetrics));
    }

    private void selectCells(final File inputBAM, final File cellListFile, final File metricsFile) {
        SelectCellsByNumTranscripts cellSelector = new SelectCellsByNumTranscripts();
        cellSelector.INPUT = inputBAM;
        cellSelector.MIN_TRANSCRIPTS_PER_CELL = 3;
        cellSelector.READ_MQ = 10;
        cellSelector.LOCUS_FUNCTION_LIST = Arrays.asList(LocusFunction.CODING, LocusFunction.UTR);
        cellSelector.OUTPUT = cellListFile;
        cellSelector.METRICS = metricsFile;

        Assert.assertEquals(cellSelector.doWork(), 0);
    }
}
