/*
 * MIT License
 *
 * Copyright 2018 Broad Institute
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
package org.broadinstitute.dropseqrna.readtrimming;

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;

public class PolyATrimmerTest {
    private static final File TESTDATA_DIR = new File("testdata/org/broadinstitute/dropseq/readtrimming");
    private static final File INPUT = new File(TESTDATA_DIR, "N701.subset.tagged_filtered_start_seq_trimmed.sam");

    // There are already tests of the real work, so just confirm that CLP runs to completion.
    @Test(dataProvider = "testClpDataProvider")
    public void testClp(final boolean newTrimmer) throws IOException {
        final File tempDir = Files.createTempDirectory("PolyATrimmerTest.").toFile();
        final Log.LogLevel saveLogLevel = Log.getGlobalLogLevel();
        Log.setGlobalLogLevel(Log.LogLevel.DEBUG);
        try {
            final PolyATrimmer clp = new PolyATrimmer();
            clp.INPUT = INPUT;
            clp.OUTPUT = File.createTempFile("PolyATrimmerTest.", ".sam");
            clp.OUTPUT.deleteOnExit();
            clp.OUTPUT_SUMMARY = File.createTempFile("PolyATrimmerTest.", ".summary");
            clp.OUTPUT_SUMMARY.deleteOnExit();
            clp.TMP_DIR = Collections.singletonList(tempDir);
            tempDir.deleteOnExit();
            clp.MISMATCHES = 0;
            clp.NUM_BASES = 6;
            clp.VALIDATION_STRINGENCY = ValidationStringency.STRICT;
            clp.USE_NEW_TRIMMER = newTrimmer;
            Assert.assertEquals(clp.doWork(), 0);
            final File expectedResult = new File(TESTDATA_DIR, String.format("N701.%s_trimmer.sam", newTrimmer? "new": "old"));
            TestUtils.assertSamFilesSame(clp.OUTPUT, expectedResult);
        } finally {
            Log.setGlobalLogLevel(saveLogLevel);
            IOUtil.recursiveDelete(tempDir.toPath());
        }
    }

    @DataProvider(name="testClpDataProvider")
    public Object[][]testClpDataProvider() {
        return new Object[][]{
                new Object[]{true},
                new Object[]{false},
        };
    }

}
