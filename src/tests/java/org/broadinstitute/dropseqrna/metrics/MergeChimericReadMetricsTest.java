/*
 * MIT License
 *
 * Copyright 2022 Broad Institute
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

import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.dropseqrna.barnyard.MarkChimericReads;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class MergeChimericReadMetricsTest {
 private static final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/metrics/MergeChimericReadMetrics");
 private static final List<File> TEST_FILES =
         IntStream.rangeClosed(0, 2).mapToObj(i -> new File(TEST_DATA_DIR, String.format("N701.%d.chimeric_read_metrics", i))).
                 collect(Collectors.toList());

 @Test
 public void testBasic() {
  final MergeChimericReadMetrics clp = new MergeChimericReadMetrics();
  clp.INPUT = TEST_FILES;
  clp.OUTPUT = TestUtils.getTempReportFile("MergeChimericReadMetricsTest.", ".chimeric_read_metrics");
  clp.DELETE_INPUTS = false;
  Assert.assertEquals(clp.doWork(), 0);
  final List<MarkChimericReads.MarkChimericReadMetrics> outputBeans = MetricsFile.readBeans(clp.OUTPUT);
  Assert.assertEquals(outputBeans.size(), 1);
  long expectNumUmis = 0L;
  for (final File f : TEST_FILES) {
   List<MarkChimericReads.MarkChimericReadMetrics> beans = MetricsFile.readBeans(f);
   Assert.assertEquals(beans.size(), 1);
   expectNumUmis += beans.getFirst().NUM_UMIS;
  }
  Assert.assertEquals(outputBeans.getFirst().NUM_UMIS, expectNumUmis);
 }
}
