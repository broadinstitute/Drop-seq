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


import com.google.cloud.storage.contrib.nio.CloudStorageFileSystemProvider;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.nio.PicardHtsPath;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.net.URI;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class FileListParsingUtilsTest {

    @Test
    public void testFileListParsing() throws IOException {
        // In the temp dir, create 2 files matching the pattern "testFileListParsing.pattern.*.txt"
        final String PATTERN = "testFileListParsing.pattern.*.txt";
        final File patternFile1 = TestUtils.getTempReportFile("testFileListParsing.pattern.", ".txt");
        final File patternFile2 = TestUtils.getTempReportFile("testFileListParsing.pattern.", ".txt");
        final File noPatternFile1 = TestUtils.getTempReportFile("testFileListParsing.no.pattern.", ".txt");
        final File noPatternFile2 = TestUtils.getTempReportFile("testFileListParsing.no.pattern.", ".txt");

        final File tempDir = patternFile1.getParentFile();

        final File fileList = TestUtils.getTempReportFile("testFileListParsing.", ".file_list");
        final PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(fileList));
        out.println(noPatternFile2.getAbsolutePath());
        out.close();

        List<File> expandedList = FileListParsingUtils.expandFileList(Arrays.asList(new File(tempDir, PATTERN), noPatternFile1, fileList));
        final List<File> expected = Arrays.asList(
                patternFile1,
                patternFile2,
                noPatternFile1,
                noPatternFile2
        );

        Collections.sort(expected);
        Collections.sort(expandedList);
        Assert.assertEquals(expected, expandedList);
    }

    @Test
    public void testPathListParsing() throws IOException {
        // In the temp dir, create 2 files matching the pattern "testPathListParsing.pattern.*.txt"
        final String PATTERN = "testPathListParsing.pattern.*.txt";
        final Path patternFile1 =
                TestUtils.getTempReportFile("testPathListParsing.pattern.", ".txt").toPath();
        final Path patternFile2 =
                TestUtils.getTempReportFile("testPathListParsing.pattern.", ".txt").toPath();
        final Path noPatternFile1 =
                TestUtils.getTempReportFile("testPathListParsing.no.pattern.", ".txt").toPath();
        final Path noPatternFile2 =
                TestUtils.getTempReportFile("testPathListParsing.no.pattern.", ".txt").toPath();

        final Path tempDir = patternFile1.getParent();

        final Path fileList = TestUtils.getTempReportFile("testPathListParsing.", ".file_list").toPath();
        final PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(fileList));
        out.println(FileUtils.toAbsoluteString(noPatternFile2));
        out.close();

        final List<Path> expandedList =
                FileListParsingUtils.expandPathList(Arrays.asList(tempDir.resolve(PATTERN), noPatternFile1, fileList));
        final List<Path> expected = Arrays.asList(
                patternFile1,
                patternFile2,
                noPatternFile1,
                noPatternFile2
        );

        Collections.sort(expected);
        Collections.sort(expandedList);
        Assert.assertEquals(expected, expandedList);
    }

    @DataProvider(name = "toAbsoluteStringData")
    public Object[][] toAbsoluteStringData() {
        return new Object[][]{
                {null, null},
                {Paths.get("/tmp/foo"), "/tmp/foo"},
                {Paths.get(URI.create("file:/tmp/foo")), "/tmp/foo"},
                {new PicardHtsPath("/tmp/foo").toPath(), "/tmp/foo"},
                {new CloudStorageFileSystemProvider().getPath("gs://bucket/foo"), "gs://bucket/foo"},
        };
    }

    @Test(dataProvider = "toAbsoluteStringData")
    public void testToAbsoluteString(final Path path, final String expected) {
        Assert.assertEquals(FileUtils.toAbsoluteString(path), expected);
    }

    @DataProvider(name = "toPrettyStringData")
    public Object[][] toPrettyStringData() {
        return new Object[][]{
                {null, null},
                {Paths.get("/tmp/foo"), "/tmp/foo"},
                {Paths.get(URI.create("file:/tmp/foo")), "/tmp/foo"},
                {new PicardHtsPath("/tmp/foo").toPath(), "/tmp/foo"},
                {new CloudStorageFileSystemProvider().getPath("gs://bucket/foo"), "gs://bucket/foo"},
        };
    }

    @Test(dataProvider = "toPrettyStringData")
    public void testToPrettyString(final Path path, final String expected) {
        Assert.assertEquals(FileUtils.toPrettyString(path), expected);
    }
}
