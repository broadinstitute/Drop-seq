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
package org.broadinstitute.dropseqrna.utils;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.text.DecimalFormat;
import java.util.Collections;
import java.util.Map;
import java.util.function.Function;

public class TestUtilsTest
{
 private static final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/utils/TestUtilsTest");
 private static final File BASELINE = new File(TEST_DATA_DIR, "AssignCellsToSamplesVerbose.tsv");
 private static final File SIMILAR = new File(TEST_DATA_DIR, "AssignCellsToSamplesVerboseSimilar.tsv");
 private static final File DIFFERENT = new File(TEST_DATA_DIR, "AssignCellsToSamplesVerboseDifferent.tsv");
 private static final File UNTRANSFORMED_DIFFERENT = new File(TEST_DATA_DIR, "AssignCellsToSamplesVerboseUntransformedDifferent.tsv");


 @Test
 public void testIdenticalTsv() {
  Assert.assertTrue(TestUtils.testFilesSame(BASELINE, BASELINE));
 }

 final DecimalFormat Format12DecimalPlaces = new DecimalFormat("0.############");

 final Map<Integer, Function<String, Object>> VerboseOutputTransformerMap = Collections.singletonMap(
         16, s -> Format12DecimalPlaces.format(Double.valueOf(s))
 );
 @Test
 public void testCloseEnoughTsv() {
  Assert.assertFalse(TestUtils.testFilesSame(BASELINE, SIMILAR));
  TestUtils.testTabularFilesSame(BASELINE, SIMILAR, VerboseOutputTransformerMap);
 }

 @Test()
 public void testDifferentTsv() {
  Assert.assertFalse(TestUtils.testTabularFilesSame(BASELINE, DIFFERENT, VerboseOutputTransformerMap));
 }
 @Test()
 public void testUntransformedDifferentTsv() {
  Assert.assertFalse(TestUtils.testFilesSame(BASELINE, UNTRANSFORMED_DIFFERENT));
  Assert.assertFalse(TestUtils.testTabularFilesSame(BASELINE, UNTRANSFORMED_DIFFERENT, VerboseOutputTransformerMap));
 }
}
