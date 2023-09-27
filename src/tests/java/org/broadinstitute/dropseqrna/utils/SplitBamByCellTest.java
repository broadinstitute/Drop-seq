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

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.testng.Assert;
import org.testng.Reporter;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.analysis.CollectAlignmentSummaryMetrics;
import picard.sam.MergeSamFiles;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.*;
import java.util.stream.Collectors;


public class SplitBamByCellTest {

    public static final File TEST_BAM = new File("testdata/org/broadinstitute/dropseq/barnyard/DgeStrandFuncTest/DgeStrandFuncTest.bam");
    private static final File EXPECTED_REPORT = new File ("testdata/org/broadinstitute/dropseq/utils/SplitBamByCell.report");
    private static final File EXPECTED_MANIFEST = new File ("testdata/org/broadinstitute/dropseq/utils/SplitBamByCell.split_bam_manifest");
    private static final File QUERYNAME_SORTED_BAM = new File("testdata/org/broadinstitute/dropseq/metrics/compute_umi_sharing.unmapped.sam");
    private static final String GENE_FUNCTION_TAG = "gf";

    @Test
	public void testDoWork() {
        final SplitBamByCell bamSplitter = new SplitBamByCell();

        File outputBAMList = TestUtils.getTempReportFile("SplitBamByCell.", ".list");
        File mergedOutputBAM = TestUtils.getTempReportFile("SplitBamByCell.", ".bam");
        File report = TestUtils.getTempReportFile("SplitBamByCell.", ".report");
        File originalMetrics = TestUtils.getTempReportFile("SplitBamByCell.", ".metrics");
        File mergedMetrics = TestUtils.getTempReportFile("SplitBamByCell.", ".metrics");

        bamSplitter.INPUT= Collections.singletonList(TEST_BAM);
        bamSplitter.OUTPUT=TestUtils.getTempReportFile("SplitBamByCell.", "." + bamSplitter.OUTPUT_SLUG + ".bam");
        bamSplitter.OUTPUT_LIST=outputBAMList;
        bamSplitter.NUM_OUTPUTS = 3;
        bamSplitter.DELETE_INPUTS = false;
        bamSplitter.OVERWRITE_EXISTING = true;
        bamSplitter.REPORT = report;
        bamSplitter.OUTPUT_MANIFEST = TestUtils.getTempReportFile("SplitBamByCell.", ".split_bam_manifest");
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
            Assert.assertEquals(fileMerger.instanceMain(args.toArray(new String[0])), 0);

            // Metrics for the input test BAM file
            args = new ArrayList<>();
            args.add("INPUT=" + TEST_BAM.getAbsolutePath());
            args.add("OUTPUT=" + originalMetrics.getAbsolutePath());
            args.add("USE_JDK_DEFLATER=" + TestUtils.isMacOs());
            final CollectAlignmentSummaryMetrics originalMetricsCollector = new CollectAlignmentSummaryMetrics();
            Assert.assertEquals(originalMetricsCollector.instanceMain(args.toArray(new String[0])), 0);

            // Metrics for the split and merged test BAM file
            args = new ArrayList<>();
            args.add("INPUT=" + mergedOutputBAM.getAbsolutePath());
            args.add("OUTPUT=" + mergedMetrics.getAbsolutePath());
            args.add("USE_JDK_DEFLATER=" + TestUtils.isMacOs());
            final CollectAlignmentSummaryMetrics mergedMetricsCollector = new CollectAlignmentSummaryMetrics();
            Assert.assertEquals(mergedMetricsCollector.instanceMain(args.toArray(new String[0])), 0);

            // Make sure the input test BAM and the merged BAM files' metrics are identical
            TestUtils.testMetricsFilesEqual(originalMetrics, mergedMetrics);

            // Compare the XC and XM tags for the input test BAM and the merged BAM files
            final CompareBAMTagValues tagsComparator = new CompareBAMTagValues();
            tagsComparator.INPUT_1 = TEST_BAM;
            tagsComparator.INPUT_2 = mergedOutputBAM;
            tagsComparator.TAGS = new ArrayList<>(Arrays.asList("XC", "XM"));
            Assert.assertEquals(tagsComparator.doWork(), 0);

            Assert.assertTrue(TestUtils.testMetricsFilesEqual(report, EXPECTED_REPORT));

            // Manifest content will not be identical because of temp file names, so just compare fields expected
            // to be consistent.
            final Map<String, Integer> expectedManifest = loadCellBarcodeManifest(EXPECTED_MANIFEST);
            final Map<String, Integer> actualManifest = loadCellBarcodeManifest(bamSplitter.OUTPUT_MANIFEST);
            Assert.assertEquals(expectedManifest, actualManifest);
		} catch (IOException ex) {
			throw new RuntimeException("Error running the test", ex);
		}
    }

    private Map<String, Integer> loadCellBarcodeManifest(final File manifest) {
        final List<CellBarcodeSplitBamMetric> beans = MetricsFile.readBeans(manifest);
        return(beans.stream().collect(Collectors.toMap(b -> b.CELL_BARCODE, b -> b.SPLIT_BAM_INDEX)));
    }

    @Test
    public void testTargetSplitSizeAndOverwrite() throws IOException {
        SplitBamByCell bamSplitter;

        File inputBam = TestUtils.getTempReportFile("SplitBamByCell.", ".bam");
        File outputBAMList = TestUtils.getTempReportFile("SplitBamByCell.", ".list");
        File reportFile = TestUtils.getTempReportFile("SplitBamByCell.", ".split_bam_report");

        Files.copy(TEST_BAM.toPath(), inputBam.toPath(), StandardCopyOption.REPLACE_EXISTING);

        bamSplitter = new SplitBamByCell();
        bamSplitter.INPUT = Collections.singletonList(inputBam);
        bamSplitter.OUTPUT = TestUtils.getTempReportFile("SplitBamByCell.", "." + bamSplitter.OUTPUT_SLUG + ".bam");
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
        final File previousOutputTemplate = bamSplitter.OUTPUT;

        // It should fail since OVERWRITE_EXISTING is false
        boolean exceptionThrown = false;
        try {
            bamSplitter = new SplitBamByCell();
            bamSplitter.INPUT = Collections.singletonList(inputBam);
            bamSplitter.OUTPUT = previousOutputTemplate;
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
        bamSplitter.INPUT = Collections.singletonList(inputBam);
        bamSplitter.OUTPUT = previousOutputTemplate;
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

    @Test
    public void testSplitEqually() {
        final SplitBamByCell bamSplitter = new SplitBamByCell();

        File mergedOutputBAM = TestUtils.getTempReportFile("SplitBamByCell.", ".sam");

        bamSplitter.INPUT= Collections.singletonList(QUERYNAME_SORTED_BAM);
        bamSplitter.OUTPUT=TestUtils.getTempReportFile("SplitBamByCell.", "." + bamSplitter.OUTPUT_SLUG + ".bam");
        bamSplitter.OUTPUT_LIST= TestUtils.getTempReportFile("SplitBamByCell.", ".list");
        bamSplitter.NUM_OUTPUTS = 3;
        bamSplitter.DELETE_INPUTS = false;
        bamSplitter.OVERWRITE_EXISTING = true;
        bamSplitter.REPORT = TestUtils.getTempReportFile("SplitBamByCell.", ".report");
        bamSplitter.SPLIT_TAG = null;
        TestUtils.setInflaterDeflaterIfMacOs();
        Assert.assertEquals(bamSplitter.doWork(), 0);

        List<File> splitBAMFileList = FileListParsingUtils.readFileList(bamSplitter.OUTPUT_LIST);
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
        args.add("SORT_ORDER=queryname");
        Assert.assertEquals(fileMerger.instanceMain(args.toArray(new String[0])), 0);
        TestUtils.assertSamRecordsSame(mergedOutputBAM, QUERYNAME_SORTED_BAM);
    }

    @Test
    public void testSplitGeneThresholdDefault() {
        final SplitBamByCell bamSplitter = new SplitBamByCell();

        final File mergedOutputBAM = TestUtils.getTempReportFile("SplitBamByCell.merged", ".sam");
        final File filterTestBAM = TestUtils.getTempReportFile("SplitBamByCell.filtered", ".sam");
        mergedOutputBAM.deleteOnExit();
        filterTestBAM.deleteOnExit();

        bamSplitter.INPUT = Collections.singletonList(TEST_BAM);
        bamSplitter.OUTPUT =
                TestUtils.getTempReportFile("SplitBamByCell.", "." + bamSplitter.OUTPUT_SLUG + ".bam");
        bamSplitter.OUTPUT_LIST = TestUtils.getTempReportFile("SplitBamByCell.", ".list");
        bamSplitter.NUM_OUTPUTS = 3;
        bamSplitter.DELETE_INPUTS = false;
        bamSplitter.OVERWRITE_EXISTING = true;
        bamSplitter.REPORT = TestUtils.getTempReportFile("SplitBamByCell.", ".report");
        bamSplitter.SPLIT_TAG = GENE_FUNCTION_TAG;
        TestUtils.setInflaterDeflaterIfMacOs();
        Assert.assertEquals(bamSplitter.doWork(), 0);

        final List<File> splitBAMFileList = FileListParsingUtils.readFileList(bamSplitter.OUTPUT_LIST);
        for (final File f : splitBAMFileList) {
            f.deleteOnExit();
        }

        // Merge the split BAM files
        final MergeSamFiles fileMerger = new MergeSamFiles();
        final List<String> fileMergerArgs = new ArrayList<>();
        for (final File splitBAMFilePath : splitBAMFileList) {
            fileMergerArgs.add("INPUT=" + splitBAMFilePath.getAbsolutePath());
        }
        fileMergerArgs.add("OUTPUT=" + mergedOutputBAM.getAbsolutePath());
        fileMergerArgs.add("USE_JDK_DEFLATER=" + TestUtils.isMacOs());
        fileMergerArgs.add("SORT_ORDER=coordinate");

        Assert.assertEquals(fileMerger.instanceMain(fileMergerArgs.toArray(new String[0])), 0);

        // Filter the original BAM file (Uses a little bit of code+cpu here versus checking in megabytes of data to git)
        final FilterBam fileFilter = new FilterBam();
        final List<String> fileFilterArgs = new ArrayList<>();
        fileFilterArgs.add("INPUT=" + TEST_BAM.getAbsolutePath());
        fileFilterArgs.add("OUTPUT=" + filterTestBAM.getAbsolutePath());
        fileFilterArgs.add("USE_JDK_DEFLATER=" + TestUtils.isMacOs());
        fileFilterArgs.add("TAG_RETAIN=" + GENE_FUNCTION_TAG);

        Assert.assertEquals(fileFilter.instanceMain(fileFilterArgs.toArray(new String[0])), 0);

        TestUtils.assertSamRecordsSame(mergedOutputBAM, filterTestBAM);
    }

    @Test
    public void testSplitGeneThresholdBadTag() {
        final SplitBamByCell bamSplitter = new SplitBamByCell();

        final File mergedOutputBAM = TestUtils.getTempReportFile("SplitBamByCell.merged.", ".sam");
        mergedOutputBAM.deleteOnExit();

        bamSplitter.INPUT = Collections.singletonList(TEST_BAM);
        bamSplitter.OUTPUT =
                TestUtils.getTempReportFile("SplitBamByCell.", "." + bamSplitter.OUTPUT_SLUG + ".bam");
        bamSplitter.OUTPUT_LIST = TestUtils.getTempReportFile("SplitBamByCell.", ".list");
        bamSplitter.NUM_OUTPUTS = 3;
        bamSplitter.DELETE_INPUTS = false;
        bamSplitter.OVERWRITE_EXISTING = true;
        bamSplitter.REPORT = TestUtils.getTempReportFile("SplitBamByCell.", ".report");
        bamSplitter.SPLIT_TAG = "BD"; // any tag that isn't present in the test BAM
        TestUtils.setInflaterDeflaterIfMacOs();
        boolean exceptionThrown = false;
        try {
            bamSplitter.doWork();
        } catch (final IllegalArgumentException e) {
            exceptionThrown = true;
        }
        Assert.assertTrue(exceptionThrown);
    }

    @Test(dataProvider = "testSplitGeneThresholdFailDataProvider")
    public void testSplitGeneThresholdFail(final Double passingReadThreshold) {
        final SplitBamByCell bamSplitter = new SplitBamByCell();

        final File mergedOutputBAM = TestUtils.getTempReportFile("SplitBamByCell.merged.", ".sam");
        mergedOutputBAM.deleteOnExit();

        bamSplitter.INPUT = Collections.singletonList(TEST_BAM);
        bamSplitter.OUTPUT =
                TestUtils.getTempReportFile("SplitBamByCell.", "." + bamSplitter.OUTPUT_SLUG + ".bam");
        bamSplitter.OUTPUT_LIST = TestUtils.getTempReportFile("SplitBamByCell.", ".list");
        bamSplitter.NUM_OUTPUTS = 3;
        bamSplitter.DELETE_INPUTS = false;
        bamSplitter.OVERWRITE_EXISTING = true;
        bamSplitter.REPORT = TestUtils.getTempReportFile("SplitBamByCell.", ".report");
        bamSplitter.SPLIT_TAG = GENE_FUNCTION_TAG;
        bamSplitter.PASSING_READ_THRESHOLD = passingReadThreshold;
        TestUtils.setInflaterDeflaterIfMacOs();
        boolean exceptionThrown = false;
        try {
            bamSplitter.doWork();
        } catch (final IllegalArgumentException e) {
            exceptionThrown = true;
        }
        Assert.assertTrue(exceptionThrown);
    }

    @DataProvider
    public Object[][] testSplitGeneThresholdFailDataProvider() {
        return new Object[][] {
                {-1d},
                {0.9},
                {60000d},
        };
    }

    @Test(dataProvider = "testDeleteInputsAndDefaultOutputsDataProvider")
    public void testDeleteInputsAndDefaultOutputs(final Boolean DELETE_INPUT_INDICES, final boolean createIndex) throws IOException {
        final File tempDir = Files.createTempDirectory("SplitBamByCellTest.").toFile();
        tempDir.deleteOnExit();
        final File inputFile = new File(tempDir, TEST_BAM.getName());
        inputFile.deleteOnExit(); // Shouldn't be necessary, because DELETE_INPUTS==true
        IOUtil.copyFile(TEST_BAM, inputFile);
        final File index;
        if (createIndex) {
            index = createBamIndexPath(inputFile);
            Assert.assertTrue(index.createNewFile());
            index.deleteOnExit();
        } else {
            index = null;
        }
        final File symlink1 = new File(tempDir, "symlink1" + AbstractSplitBamClp.BAM_EXTENSION);
        Files.createSymbolicLink(symlink1.toPath(), inputFile.toPath());
        symlink1.deleteOnExit();
        final File indexSymlink1;
        if (createIndex) {
            indexSymlink1 = createBamIndexPath(symlink1);
            Files.createSymbolicLink(indexSymlink1.toPath(), index.toPath());
            indexSymlink1.deleteOnExit();
        } else {
            indexSymlink1 = null;
        }
        final String symlinkBasename = "symlink";
        final File inputSymlink = new File(tempDir, symlinkBasename + AbstractSplitBamClp.BAM_EXTENSION);
        Files.createSymbolicLink(inputSymlink.toPath(), symlink1.toPath());
        inputSymlink.deleteOnExit();
        if (createIndex) {
            final File indexSymlink = createBamIndexPath(inputSymlink);
            Files.createSymbolicLink(indexSymlink.toPath(), indexSymlink1.toPath());
            indexSymlink.deleteOnExit();
        }

        final SplitBamByCell bamSplitter = new SplitBamByCell();
        bamSplitter.INPUT= Collections.singletonList(inputSymlink);
        final int NUM_OUTPUTS = 3;
        bamSplitter.NUM_OUTPUTS = NUM_OUTPUTS;
        bamSplitter.DELETE_INPUTS = true;
        bamSplitter.DELETE_INPUT_INDICES = DELETE_INPUT_INDICES;
        TestUtils.setInflaterDeflaterIfMacOs();
        Assert.assertEquals(bamSplitter.doWork(), 0);
        List<File> outputFiles = Arrays.asList(tempDir.listFiles((dir, name) -> !name.startsWith(".nfs")));
        outputFiles.forEach(File::deleteOnExit);
        if (createIndex && DELETE_INPUT_INDICES != null && !DELETE_INPUT_INDICES) {
            outputFiles = outputFiles.stream().filter(f -> !f.getName().endsWith(".bai")).collect(Collectors.toList());
        }
        final Set<String> outputFileNames = outputFiles.stream().map(File::getName).collect(Collectors.toSet());
        Reporter.log(String.format("Files found: %s", StringUtil.join(", ", outputFileNames)), true);

        final List<String> expectedExtensions = new ArrayList<>(Arrays.asList(AbstractSplitBamClp.BAM_LIST_EXTENSION, AbstractSplitBamClp.BAM_REPORT_EXTENSION));
        for (int i = 0; i < NUM_OUTPUTS; ++i) {
            expectedExtensions.add("." + i + AbstractSplitBamClp.BAM_EXTENSION);
        }
        // split BAMs + bam_list and report
        Assert.assertEquals(expectedExtensions.size(), outputFiles.size());
        for (final String expectedExtension: expectedExtensions) {
            Assert.assertTrue(outputFileNames.contains(symlinkBasename + expectedExtension), "Test presence of " + expectedExtension
                    );
        }
    }

    private static File createBamIndexPath(final File bamPath) {
        String basename = bamPath.getName();
        basename = basename.substring(0, basename.length() - 1);
        return new File(bamPath.getParentFile(), basename + "i");
    }

    @DataProvider(name = "testDeleteInputsAndDefaultOutputsDataProvider")
    public Object[][] testDeleteInputsAndDefaultOutputsDataProvider() {
        return new Object[][] {
                {true, true},
                {true, false},
                {false, true},
                {false, false},
                {null, true},
                {null, false},
        };
    }
}
