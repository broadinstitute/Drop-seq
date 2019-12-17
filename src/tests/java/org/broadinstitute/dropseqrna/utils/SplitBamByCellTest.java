/*
 * MIT License
 *
 * Copyright 2019 Broad Institute
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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.junit.Assert;
import org.testng.annotations.Test;
import picard.analysis.CollectAlignmentSummaryMetrics;
import picard.sam.MergeSamFiles;


public class SplitBamByCellTest {

    private static final File TEST_BAM = new File("testdata/org/broadinstitute/dropseq/barnyard/DgeStrandFuncTest/DgeStrandFuncTest.bam");
    private static File EXPECTED_REPORT = new File ("testdata/org/broadinstitute/dropseq/utils/SplitBamByCell.report");

    @Test
	public void testDoWork() {
        final SplitBamByCell bamSplitter = new SplitBamByCell();

		File outputBAM = getTempReportFile("SplitBamByCell", bamSplitter.OUTPUT_SLUG + ".bam");
        File outputBAMList = getTempReportFile("SplitBamByCell", ".list");
        File mergedOutputBAM = getTempReportFile("SplitBamByCell", ".bam");
        File report = getTempReportFile("SplitBamByCell", ".report");
        File originalMetrics = getTempReportFile("SplitBamByCell", ".metrics");
        File mergedMetrics = getTempReportFile("SplitBamByCell", ".metrics");

        bamSplitter.INPUT=Arrays.asList(TEST_BAM);
        bamSplitter.OUTPUT=outputBAM;
        bamSplitter.OUTPUT_LIST=outputBAMList;
        bamSplitter.NUM_OUTPUTS = 3;
        bamSplitter.REPORT = report;
		Assert.assertEquals(bamSplitter.doWork(), 0);

        List<String> splitBAMFileList;
		try {
            splitBAMFileList = FileUtils.readLines(outputBAMList);
            for (String filePath : splitBAMFileList) {
                new File(filePath).deleteOnExit();
            }

            // Merge the split BAM files
            final MergeSamFiles fileMerger = new MergeSamFiles();
            List<String> args = new ArrayList<>();
            for (String splitBAMFilePath : splitBAMFileList)  {
                args.add("INPUT=" + splitBAMFilePath);
            }
            args.add("OUTPUT=" + mergedOutputBAM.getAbsolutePath());
            Assert.assertEquals(fileMerger.instanceMain(args.toArray(new String[args.size()])), 0);

            // Metrics for the input test BAM file
            args = new ArrayList<String>();
            args.add("INPUT=" + TEST_BAM.getAbsolutePath());
            args.add("OUTPUT=" + originalMetrics.getAbsolutePath());
            final CollectAlignmentSummaryMetrics originalMetricsCollector = new CollectAlignmentSummaryMetrics();
            Assert.assertEquals(originalMetricsCollector.instanceMain(args.toArray(new String[args.size()])), 0);

            // Metrics for the split and merged test BAM file
            args = new ArrayList<String>();
            args.add("INPUT=" + mergedOutputBAM.getAbsolutePath());
            args.add("OUTPUT=" + mergedMetrics.getAbsolutePath());
            final CollectAlignmentSummaryMetrics mergedMetricsCollector = new CollectAlignmentSummaryMetrics();
            Assert.assertEquals(mergedMetricsCollector.instanceMain(args.toArray(new String[args.size()])), 0);

            // Make sure the input test BAM and the merged BAM files' metrics are identical
            TestUtils.testMetricsFilesEqual(originalMetrics, mergedMetrics);

            // Compare the XC and XM tags for the input test BAM and the merged BAM files
            final CompareBAMTagValues tagsComparator = new CompareBAMTagValues();
            tagsComparator.INPUT_1 = TEST_BAM;
            tagsComparator.INPUT_2 = mergedOutputBAM;
            tagsComparator.TAGS = new ArrayList<>(Arrays.asList("XC", "XM"));
            Assert.assertEquals(tagsComparator.doWork(), 0);

			boolean t1 = FileUtils.contentEquals(report, EXPECTED_REPORT);
			Assert.assertTrue(t1);
		} catch (IOException e) {
			e.printStackTrace();
		}
    }

    private File getTempReportFile (final String prefix, final String suffix) {
        File tempFile=null;

        try {
            tempFile = File.createTempFile(prefix, suffix);
            tempFile.deleteOnExit();
        } catch (IOException e1) {
            e1.printStackTrace();
        }
        return tempFile;
    }
}
