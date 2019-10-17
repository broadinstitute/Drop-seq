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
package org.broadinstitute.dropseqrna.metrics;

import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

public class CountUnmatchedSampleIndicesTest {
    private static final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/metrics/CountUnmatchedSampleIndices");
    private static final File EXPECTED_METRICS = new File(TEST_DATA_DIR, "expected.unmatched_index_metrics");
    private static final int NUM_BARCODE_FILES = 9;
    private static final String BARCODE_FILE_TEMPLATE = "s_1_1110%d_barcode.txt.gz";

    @Test(dataProvider = "testBasicDataProvider")
    public void testBasic(int numThreads) throws IOException {
        final CountUnmatchedSampleIndices clp = new CountUnmatchedSampleIndices();
        clp.NUM_THREADS=numThreads;
        clp.OUTPUT = File.createTempFile("CountUnmatchedSampleIndicesTest.", "unmatched_index_metrics");
        clp.OUTPUT.deleteOnExit();
        clp.BARCODE_FILES = new ArrayList<>(NUM_BARCODE_FILES);
        for (int i = 1; i <= NUM_BARCODE_FILES; ++i) {
            clp.BARCODE_FILES.add(new File(TEST_DATA_DIR, String.format(BARCODE_FILE_TEMPLATE, i)));
        }
        Assert.assertEquals(clp.doWork(), 0);
        Assert.assertTrue(TestUtils.testMetricsFilesEqual(EXPECTED_METRICS, clp.OUTPUT));
    }

    @Test(dataProvider = "testBasicDataProvider")
    public void testDirectoryInput(int numThreads) throws IOException {
        final CountUnmatchedSampleIndices clp = new CountUnmatchedSampleIndices();
        clp.NUM_THREADS=numThreads;
        clp.OUTPUT = File.createTempFile("CountUnmatchedSampleIndicesTest.", "unmatched_index_metrics");
        clp.OUTPUT.deleteOnExit();
        clp.BARCODE_FILES = Arrays.asList(TEST_DATA_DIR);
        Assert.assertEquals(clp.doWork(), 0);
        Assert.assertTrue(TestUtils.testMetricsFilesEqual(EXPECTED_METRICS, clp.OUTPUT));
    }
    @DataProvider(name="testBasicDataProvider")
    public Object[][]testBasicDataProvider() {
        return new Object[][] {
            {1}, {2}, {4}, {8}
        };
    }
}
