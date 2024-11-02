/*
 * MIT License
 *
 * Copyright 2024 Broad Institute
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

import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class MergeSplitDgesTest {
 private static final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/barnyard/digitalexpression/tools");
 private static final File INPUT1 = new File(TEST_DATA_DIR, "MergeSplitDge.1.digital_expression.txt.gz");
 private static final File INPUT2 = new File(TEST_DATA_DIR, "MergeSplitDge.2.digital_expression.txt.gz");
 private static final List<File> INPUTS = Arrays.asList(INPUT1, INPUT2);
 private static final File EXPECTED_OUTPUT = new File(TEST_DATA_DIR, "MergeSplitDge.expected.digital_expression.txt.gz");
 private static final File EXPECTED_PREFIX_OUTPUT = new File(TEST_DATA_DIR, "MergeSplitDge.prefix.expected.digital_expression.txt.gz");


 @Test
 public void testBasic() {
  final MergeSplitDges clp = new MergeSplitDges();
  clp.INPUT = INPUTS;
  clp.OUTPUT = TestUtils.getTempReportFile("MergeDge", ".txt.gz");
  Assert.assertEquals(clp.doWork(), 0);
  Assert.assertTrue(TestUtils.dgeMatricesAreEqual(clp.OUTPUT, EXPECTED_OUTPUT));
 }

 @Test(expectedExceptions = IllegalArgumentException.class)
 public void testBarcodeCollision() {
  final MergeSplitDges clp = new MergeSplitDges();
  clp.INPUT = Arrays.asList(INPUT1, INPUT1);
  clp.OUTPUT = TestUtils.getTempReportFile("MergeDge", ".txt.gz");
  clp.doWork();
 }

 @Test
 void testPrefix() {
  final MergeSplitDges clp = new MergeSplitDges();
  clp.INPUT = Arrays.asList(INPUT1, INPUT1);
  clp.PREFIX = Arrays.asList("1_", "2_");
  clp.OUTPUT = TestUtils.getTempReportFile("MergeDge", ".txt.gz");
  Assert.assertEquals(clp.doWork(), 0);
  Assert.assertTrue(TestUtils.dgeMatricesAreEqual(clp.OUTPUT, EXPECTED_PREFIX_OUTPUT));
 }

 @Test
 public void testEmptyPrefix() {
  final MergeSplitDges clp = new MergeSplitDges();
  clp.INPUT = INPUTS;
  clp.PREFIX = Collections.emptyList();
  clp.OUTPUT = TestUtils.getTempReportFile("MergeDge", ".txt.gz");
  Assert.assertEquals(clp.doWork(), 0);
  Assert.assertTrue(TestUtils.dgeMatricesAreEqual(clp.OUTPUT, EXPECTED_OUTPUT));
 }

 @Test(expectedExceptions = IllegalArgumentException.class)
 void testOutOfOrderGenes() {
  final MergeSplitDges clp = new MergeSplitDges();
  clp.INPUT = List.of(new File(TEST_DATA_DIR, "MergeSplitDge.wrongOrder.digital_expression.txt.gz"));
  clp.OUTPUT = TestUtils.getTempReportFile("MergeDge", ".txt.gz");
  clp.doWork();

 }
}
