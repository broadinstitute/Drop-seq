/*
 * MIT License
 *
 * Copyright 2023 Broad Institute
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

import htsjdk.samtools.util.RuntimeIOException;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class MergeUMIReadIntervalsTest {

 private static final File SECOND_FILE_TO_MERGE = new File(GatherUMIReadIntervalsTest.TEST_DATA_DIR, "GatherUMIReadIntervalsTest.2.tsv");
 public static final List<File> INPUT_FILES = Arrays.asList(GatherUMIReadIntervalsTest.EXPECTED_REPORT, SECOND_FILE_TO_MERGE);
 private static final File EXPECTED_MERGED_REPORT = new File(GatherUMIReadIntervalsTest.TEST_DATA_DIR, "MergeUMIReadIntervalsTest.tsv");

 @Test
 public void testBasic() {
  testHelper(INPUT_FILES, false);
 }

 @Test
 public void testDeleteInputs() {
  final List<File> copiedInputFiles = INPUT_FILES.stream().map(f ->
  {
   final File copiedFile = TestUtils.getTempReportFile("MergeUMIReadIntervalsTest.", ".tsv");
   try {
    Files.copy(f.toPath(), copiedFile.toPath(), StandardCopyOption.REPLACE_EXISTING);
   } catch (IOException e) {
    throw new RuntimeIOException(e);
   }
   return copiedFile;
  }).collect(Collectors.toList());
  testHelper(copiedInputFiles, true);
  copiedInputFiles.forEach(f -> Assert.assertFalse(f.exists()));
 }

 private void testHelper(final List<File> inputs, boolean deleteInputs) {
  final MergeUMIReadIntervals clp = new MergeUMIReadIntervals();
  clp.INPUT = inputs;
  clp.OUTPUT = TestUtils.getTempReportFile("MergeUMIReadIntervalsTest.", ".tsv");
  clp.DELETE_INPUTS = deleteInputs;
  Assert.assertEquals(clp.doWork(), 0);
  Assert.assertTrue(TestUtils.testFilesSame(EXPECTED_MERGED_REPORT, clp.OUTPUT));

 }
}
