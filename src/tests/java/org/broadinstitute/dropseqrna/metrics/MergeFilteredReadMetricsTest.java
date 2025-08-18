/*
 * MIT License
 *
 * Copyright 2021 Broad Institute
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
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public class MergeFilteredReadMetricsTest {
 private static final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/metrics/MergeFilteredReadMetrics");
 private static final String[] INPUT_FILENAMES = {
         "input1.filtered_read_metrics",
         "input2.filtered_read_metrics",
         "input3.filtered_read_metrics"
 };
 private static final File EXPECTED_RESULT = new File(TEST_DATA_DIR, "merged.filtered_read_metrics");

 @Test
 public void testNoDelete() throws FileNotFoundException {
  testHelper(Arrays.stream(INPUT_FILENAMES).map(name -> new File(TEST_DATA_DIR, name)).collect(Collectors.toList()), false);
 }

 @Test
 public void testDelete() throws IOException {
  final List<File> inputCopies = Arrays.stream(INPUT_FILENAMES).map(name -> new File(TEST_DATA_DIR, name)).
          map(this::copyToTemp).collect(Collectors.toList());
  testHelper(inputCopies, true);
  inputCopies.forEach(f -> Assert.assertFalse(f.exists()));
 }

 private File copyToTemp(final File source) {
  final File target = TestUtils.getTempReportFile("MergeFilteredReadMetrics.", ".filtered_read_metrics");
  try {
   Files.copy(source.toPath(), target.toPath(), StandardCopyOption.REPLACE_EXISTING);
  } catch (IOException e) {
   throw new RuntimeIOException(e);
  }
  return target;
 }

 private void testHelper(final List<File> inputs, final boolean deleteInputs) throws FileNotFoundException {
  final MergeFilteredReadMetrics clp = new MergeFilteredReadMetrics();
  clp.INPUT = inputs;
  clp.DELETE_INPUTS = deleteInputs;
  clp.OUTPUT = TestUtils.getTempReportFile("MergeFilteredReadMetrics.", ".filtered_read_metrics");
  Assert.assertEquals(clp.doWork(), 0);
  Assert.assertTrue(TestUtils.testMetricsFilesEqual(EXPECTED_RESULT, clp.OUTPUT),
          String.format("%s != %s", EXPECTED_RESULT, clp.OUTPUT));

 }
}
