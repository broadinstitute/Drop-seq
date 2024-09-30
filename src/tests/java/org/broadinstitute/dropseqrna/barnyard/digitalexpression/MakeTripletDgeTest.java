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

import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;

public class MakeTripletDgeTest {
    public static final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/barnyard/digitalexpression/MakeTripletDgeTest");

    @Test
    public void testBasic() {
        final MakeTripletDge makeTripletDge = runIt("N701.manifest.yaml");
        testIt(makeTripletDge, "N701.barcodes.tsv.gz", "N701.features.tsv.gz", "N701.matrix.mtx.gz");
    }

    @Test
    public void testBarcodeList() {
        final MakeTripletDge makeTripletDge = runIt("N701.barcode_list.manifest.yaml");
        testIt(makeTripletDge, "N701.barcode_list.barcodes.tsv.gz", "N701.features.tsv.gz", "N701.barcode_list.matrix.mtx.gz");
    }

    @Test
    public void testNoPrefix() {
        final MakeTripletDge makeTripletDge = runIt("N701.no_prefix.manifest.yaml");
        testIt(makeTripletDge, "N701.no_prefix.barcodes.tsv.gz", "N701.features.tsv.gz", "N701.matrix.mtx.gz");
    }

    @Test
    public void testTwoInputs() {
        final MakeTripletDge makeTripletDge = runIt("two_inputs.manifest.yaml");
        testIt(makeTripletDge, "two_inputs.barcodes.tsv.gz", "two_inputs.features.tsv.gz", "two_inputs.matrix.mtx.gz");
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBarcodeCollision() {
        runIt("barcode_collision.manifest.yaml");
    }

    private MakeTripletDge runIt(final String manifestName) {
        final MakeTripletDge makeTripletDge = new MakeTripletDge();
        makeTripletDge.TMP_DIR = Arrays.asList(TestUtils.createTempDirectory("MakeTripletDgeTest."));

        makeTripletDge.MANIFEST = new File(TEST_DATA_DIR, manifestName);
        makeTripletDge.OUTPUT_CELLS = TestUtils.getTempReportFile("MakeTripletDgeTest.", ".barcodes.tsv.gz");
        makeTripletDge.OUTPUT_FEATURES = TestUtils.getTempReportFile("MakeTripletDgeTest.", ".features.tsv.gz");
        makeTripletDge.OUTPUT = TestUtils.getTempReportFile("MakeTripletDgeTest.", ".matrix.mtx.gz");
        Assert.assertEquals(makeTripletDge.doWork(), 0);
        return makeTripletDge;
    }

    private void testIt(final MakeTripletDge makeTripletDge, final String barcodeFile, final String featureFile, final String matrixFile) {
        TestUtils.testFilesSame(makeTripletDge.OUTPUT_CELLS, new File(TEST_DATA_DIR, barcodeFile));
        TestUtils.testFilesSame(makeTripletDge.OUTPUT_FEATURES, new File(TEST_DATA_DIR, featureFile));
        // This file is not tabular, so we can't use TestUtils.testFilesSame
        IOUtil.assertFilesEqual(makeTripletDge.OUTPUT, new File(TEST_DATA_DIR, matrixFile));
    }
}
