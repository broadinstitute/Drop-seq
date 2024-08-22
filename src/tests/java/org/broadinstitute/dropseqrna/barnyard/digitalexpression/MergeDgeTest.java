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
package org.broadinstitute.dropseqrna.barnyard.digitalexpression;

import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.dropseqrna.barnyard.DigitalExpression;
import org.broadinstitute.dropseqrna.barnyard.SelectCellsByNumTranscripts;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.tools.DGEMatrix;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.*;

public class MergeDgeTest {
    public static final File TEST_BAM = new File("testdata/org/broadinstitute/dropseq/barnyard/DgeStrandFuncTest/DgeStrandFuncTest.bam");
    public static final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/barnyard/digitalexpression/tools");
    public static final File TEST_DGE = new File(TEST_DATA_DIR, "N701.auto.digital_expression.txt.gz");

    @Test
    public void testDoWork() {
        File cellsFile1 = TestUtils.getTempReportFile("MergeDge", ".txt.gz");
        File cellsFile2 = TestUtils.getTempReportFile("MergeDge", ".txt.gz");
        File originalCellsFile = TestUtils.getTempReportFile("MergeDge", ".txt.gz");
        File output1 = TestUtils.getTempReportFile("MergeDge", ".txt");
        File output2 = TestUtils.getTempReportFile("MergeDge", ".txt");
        File originalOutput = TestUtils.getTempReportFile("MergeDge", ".txt");
        File mergedOutput = TestUtils.getTempReportFile("MergeDge", ".txt");
        File summary1 = TestUtils.getTempReportFile("MergeDge", ".summary");
        File summary2 = TestUtils.getTempReportFile("MergeDge", ".summary");
        File originalSummary = TestUtils.getTempReportFile("MergeDge", ".summary");
        File mergedSummary = TestUtils.getTempReportFile("MergeDge", ".summary");

        List<File> splitBAMFileList = TestUtils.splitBamFile(TEST_BAM, 2);
        Assert.assertEquals(splitBAMFileList.size(), 2);

        selectCells(TEST_BAM, originalCellsFile);
        selectCells(splitBAMFileList.get(0), cellsFile1);
        selectCells(splitBAMFileList.get(1), cellsFile2);

        // Compute DigitalExpression matrices and metrics for TEST_BAM and each of the split BAM files
        collectDigitalExpression(TEST_BAM, originalCellsFile, originalOutput, originalSummary);
        collectDigitalExpression(splitBAMFileList.get(0), cellsFile1, output1, summary1);
        collectDigitalExpression(splitBAMFileList.get(1), cellsFile2, output2, summary2);

        // Merge DigitalExpression matrices
        MergeDge dgeMatrixMerger = new MergeDge();
        dgeMatrixMerger.INPUT = Arrays.asList(output1, output2);
        dgeMatrixMerger.HEADER_STRINGENCY = DgeHeaderMerger.Stringency.LENIENT;
        dgeMatrixMerger.OUTPUT = mergedOutput;
        dgeMatrixMerger.OUTPUT_HEADER = true;
        Assert.assertEquals(dgeMatrixMerger.doWork(), 0);

        Assert.assertTrue(TestUtils.dgeMatricesAreEqual(mergedOutput, originalOutput));

        // Merge DigitalExpression summaries
        MergeDgeSummaries summaryMerger = new MergeDgeSummaries();
        summaryMerger.INPUT = Arrays.asList(summary1, summary2);
        summaryMerger.OUTPUT = mergedSummary;
        Assert.assertEquals(summaryMerger.doWork(), 0);

        Assert.assertTrue(MetricsFile.areMetricsEqual(mergedSummary, originalSummary));
    }

    private void selectCells(File inputBAM, File cellsListFile) {
        SelectCellsByNumTranscripts cellSelector = new SelectCellsByNumTranscripts();

        String[] args = new String[] {
                "INPUT=" + inputBAM.getAbsolutePath(),
                "MIN_TRANSCRIPTS_PER_CELL=1",
                "OUTPUT=" + cellsListFile.getAbsolutePath()
        };

        Assert.assertEquals(cellSelector.instanceMain(args), 0);
    }

    private void collectDigitalExpression(File inputBAM, File cellsListFile, File outputFile, File summaryFile) {
        DigitalExpression metricsCollector = new DigitalExpression();
        String[] args = new String[] {
                "INPUT=" + inputBAM.getAbsolutePath(),
                "CELL_BC_FILE=" + cellsListFile.getAbsolutePath(),
                "MIN_BC_READ_THRESHOLD=0",
                "OUTPUT_HEADER=true",
                "UNIQUE_EXPERIMENT_ID=foo",
                "OUTPUT=" + outputFile.getAbsolutePath(),
                "SUMMARY=" + summaryFile.getAbsolutePath()
        };
        Assert.assertEquals(metricsCollector.instanceMain(args), 0);
    }

    @Test
    public void testDropSeqSparse() {
        MergeDge dgeMatrixMerger = new MergeDge();
        dgeMatrixMerger.INPUT = Arrays.asList(TEST_DGE);
        dgeMatrixMerger.OUTPUT = TestUtils.getTempReportFile("MergeDge", ".DropSeqSparse.txt.gz");;
        dgeMatrixMerger.OUTPUT_HEADER = false;
        dgeMatrixMerger.OUTPUT_FORMAT = DGEMatrix.FileFormat.MM_SPARSE;
        Assert.assertEquals(dgeMatrixMerger.doWork(), 0);
        Assert.assertEquals(DGEMatrix.FileFormat.MM_SPARSE, DGEMatrix.detectFileFormat(dgeMatrixMerger.OUTPUT));

        Assert.assertTrue(TestUtils.dgeMatricesAreEqual(TEST_DGE, dgeMatrixMerger.OUTPUT));

    }

    @Test
    public void TestMMSparse10XWithGenes() {
        MergeDge dgeMatrixMerger = new MergeDge();
        dgeMatrixMerger.INPUT = Arrays.asList(TEST_DGE);
        dgeMatrixMerger.OUTPUT = TestUtils.getTempReportFile("MergeDge", ".matrix.txt.gz");;
        dgeMatrixMerger.OUTPUT_GENES = TestUtils.getTempReportFile("MergeDge", ".genes.tsv");
        dgeMatrixMerger.OUTPUT_CELLS = TestUtils.getTempReportFile("MergeDge", ".barcodes.tsv");
        dgeMatrixMerger.OUTPUT_HEADER = false;
        dgeMatrixMerger.OUTPUT_FORMAT = DGEMatrix.FileFormat.MM_SPARSE_10X;
        Assert.assertEquals(dgeMatrixMerger.doWork(), 0);
        Assert.assertEquals(DGEMatrix.detectFileFormat(dgeMatrixMerger.OUTPUT), DGEMatrix.FileFormat.MM_SPARSE_10X);

        Assert.assertTrue(TestUtils.dgeMatricesAreEqual(DGEMatrix.parseFile(TEST_DGE),
                DGEMatrix.parseFile(dgeMatrixMerger.OUTPUT, dgeMatrixMerger.OUTPUT_GENES, dgeMatrixMerger.OUTPUT_CELLS, "")));

    }

    @Test
    public void TestMMSparse10XWithFeatures() {
        MergeDge dgeMatrixMerger = new MergeDge();
        dgeMatrixMerger.INPUT = Arrays.asList(TEST_DGE);
        dgeMatrixMerger.OUTPUT = TestUtils.getTempReportFile("MergeDge", ".matrix.txt.gz");;
        dgeMatrixMerger.OUTPUT_FEATURES = TestUtils.getTempReportFile("MergeDge", ".features.tsv.gz");
        dgeMatrixMerger.OUTPUT_CELLS = TestUtils.getTempReportFile("MergeDge", ".barcodes.tsv.gz");
        dgeMatrixMerger.OUTPUT_HEADER = false;
        dgeMatrixMerger.OUTPUT_FORMAT = DGEMatrix.FileFormat.MM_SPARSE_10X;
        Assert.assertEquals(dgeMatrixMerger.doWork(), 0);
        Assert.assertEquals(DGEMatrix.detectFileFormat(dgeMatrixMerger.OUTPUT), DGEMatrix.FileFormat.MM_SPARSE_10X);

        Assert.assertTrue(TestUtils.dgeMatricesAreEqual(DGEMatrix.parseFile(TEST_DGE),
                DGEMatrix.parseFile(dgeMatrixMerger.OUTPUT, dgeMatrixMerger.OUTPUT_FEATURES, dgeMatrixMerger.OUTPUT_CELLS, "")));

    }
}
