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
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import htsjdk.samtools.util.IOUtil;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.junit.Assert;
import org.omg.PortableInterceptor.INACTIVE;
import org.testng.annotations.Test;
import picard.analysis.CollectAlignmentSummaryMetrics;
import picard.sam.MergeSamFiles;

import javax.print.DocFlavor;


public class SplitBamByCellTest {

    public static final File TEST_BAM = new File("testdata/org/broadinstitute/dropseq/barnyard/DgeStrandFuncTest/DgeStrandFuncTest.bam");
    private static File EXPECTED_REPORT = new File ("testdata/org/broadinstitute/dropseq/utils/SplitBamByCell.report");

    @Test
	public void testDoWork() {
        final SplitBamByCell bamSplitter = new SplitBamByCell();

        File outputBAMList = TestUtils.getTempReportFile("SplitBamByCell.", ".list");
        File mergedOutputBAM = TestUtils.getTempReportFile("SplitBamByCell.", ".bam");
        File report = TestUtils.getTempReportFile("SplitBamByCell.", ".report");
        File originalMetrics = TestUtils.getTempReportFile("SplitBamByCell.", ".metrics");
        File mergedMetrics = TestUtils.getTempReportFile("SplitBamByCell.", ".metrics");

        bamSplitter.INPUT=Arrays.asList(TEST_BAM);
        bamSplitter.OUTPUT=new File("SplitBamByCell." + bamSplitter.OUTPUT_SLUG + ".bam");
        bamSplitter.OUTPUT_LIST=outputBAMList;
        bamSplitter.NUM_OUTPUTS = 3;
        bamSplitter.DELETE_INPUTS = false;
        bamSplitter.OVERWRITE_EXISTING = true;
        bamSplitter.REPORT = report;
        TestUtils.setInflaterDeflaterIfMacOs();
		Assert.assertEquals(bamSplitter.doWork(), 0);

        List<File> splitBAMFileList = FileListParsingUtils.readFileList(outputBAMList);
		try {
            for (File f : splitBAMFileList) {
                f.deleteOnExit();
            }

            // Merge the split BAM files
            final MergeSamFiles fileMerger = new MergeSamFiles();
            List<String> args = new ArrayList<>();
            for (File splitBAMFilePath : splitBAMFileList)  {
                args.add("INPUT=" + splitBAMFilePath.getAbsolutePath());
            }
            args.add("OUTPUT=" + mergedOutputBAM.getAbsolutePath());
            args.add("USE_JDK_DEFLATER=" + TestUtils.isMacOs());
            Assert.assertEquals(fileMerger.instanceMain(args.toArray(new String[args.size()])), 0);

            // Metrics for the input test BAM file
            args = new ArrayList<String>();
            args.add("INPUT=" + TEST_BAM.getAbsolutePath());
            args.add("OUTPUT=" + originalMetrics.getAbsolutePath());
            args.add("USE_JDK_DEFLATER=" + TestUtils.isMacOs());
            final CollectAlignmentSummaryMetrics originalMetricsCollector = new CollectAlignmentSummaryMetrics();
            Assert.assertEquals(originalMetricsCollector.instanceMain(args.toArray(new String[args.size()])), 0);

            // Metrics for the split and merged test BAM file
            args = new ArrayList<String>();
            args.add("INPUT=" + mergedOutputBAM.getAbsolutePath());
            args.add("OUTPUT=" + mergedMetrics.getAbsolutePath());
            args.add("USE_JDK_DEFLATER=" + TestUtils.isMacOs());
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

			Assert.assertTrue(FileUtils.contentEquals(report, EXPECTED_REPORT));
		} catch (IOException ex) {
			throw new RuntimeException("Error running the test", ex);
		}
    }

    @Test
    public void testTargetSplitSizeAndOverwrite() throws IOException {
        SplitBamByCell bamSplitter;

        File inputBam = TestUtils.getTempReportFile("SplitBamByCell.", ".bam");
        File outputBAMList = TestUtils.getTempReportFile("SplitBamByCell.", ".list");
        File reportFile = TestUtils.getTempReportFile("SplitBamByCell.", ".split_bam_report");

        Files.copy(TEST_BAM.toPath(), inputBam.toPath(), StandardCopyOption.REPLACE_EXISTING);

        bamSplitter = new SplitBamByCell();
        bamSplitter.INPUT = Arrays.asList(inputBam);
        bamSplitter.OUTPUT = new File("SplitBamByCell2." + bamSplitter.OUTPUT_SLUG + ".bam");
        bamSplitter.OUTPUT_LIST = outputBAMList;
        bamSplitter.REPORT = reportFile;
        bamSplitter.TARGET_BAM_SIZE = "1.5M";
        bamSplitter.DELETE_INPUTS = false;
        Assert.assertEquals(bamSplitter.doWork(), 0);

        List<File> splitBAMFileList = FileListParsingUtils.readFileList(outputBAMList);
        Assert.assertEquals(splitBAMFileList.size(), 3);
        for (File f : splitBAMFileList) {
            f.deleteOnExit();
        }

        // It should fail since OVERWRITE_EXISTING is false
        boolean exceptionThrown = false;
        try {
            bamSplitter = new SplitBamByCell();
            bamSplitter.INPUT = Arrays.asList(inputBam);
            bamSplitter.OUTPUT = new File("SplitBamByCell." + bamSplitter.OUTPUT_SLUG + ".bam");
            bamSplitter.OUTPUT_LIST = outputBAMList;
            bamSplitter.DELETE_INPUTS = false;
            bamSplitter.TARGET_BAM_SIZE = "2M";
            bamSplitter.doWork();
        } catch (IllegalArgumentException ex) {
            exceptionThrown = true;
        }
        Assert.assertTrue(exceptionThrown);
        Assert.assertTrue(inputBam.exists());

        // Now it should overwrite the output BAM files and then delete the input BAM
        bamSplitter = new SplitBamByCell();
        bamSplitter.INPUT = Arrays.asList(inputBam);
        bamSplitter.OUTPUT = new File("SplitBamByCell." + bamSplitter.OUTPUT_SLUG + ".bam");
        bamSplitter.OUTPUT_LIST = outputBAMList;
        bamSplitter.REPORT = reportFile;
        bamSplitter.DELETE_INPUTS = false;
        bamSplitter.TARGET_BAM_SIZE = "2M";
        bamSplitter.OVERWRITE_EXISTING = true;
        bamSplitter.DELETE_INPUTS = true;
        Assert.assertEquals(bamSplitter.doWork(), 0);
        Assert.assertFalse(inputBam.exists());
    }

    @Test
    public void testReadBamList() throws IOException {
        // Create a bam_list with one relative path and one absolute, and confirm that reader
        // resolves them correctly.
        final File bamList = TestUtils.getTempReportFile("testReadBamList.", ".bam_list");
        final PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(bamList));
        final String relative = "foo.bam";
        out.println(relative);
        final String absolute = "/an/absolute/path.bam";
        out.println(absolute);
        out.close();
        final List<File> expected = Arrays.asList(
          new File(bamList.getCanonicalFile().getParent(), relative),
          new File(absolute)
        );
        Assert.assertEquals(expected, FileListParsingUtils.readFileList(bamList));
        // Confirm that when reading a symlink to a bam_list, relative paths are resolved relative to the directory
        // of the actually bam_list file, not the directory containing the symlink.
        final File otherDir = Files.createTempDirectory("testReadBamList").toFile();
        otherDir.deleteOnExit();
        final File symlink = new File(otherDir, "testReadBamList.bam_list");
        symlink.deleteOnExit();
        Files.createSymbolicLink(symlink.toPath(), bamList.toPath());
        Assert.assertEquals(expected, FileListParsingUtils.readFileList(symlink));
    }
}
