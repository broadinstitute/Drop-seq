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
package org.broadinstitute.dropseqrna.eqtl;

import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;

public class MakeMetacellsFromTripletDgeTest {
    private final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/eqtl/MakeMetacellsFromTripletDgeTest");

    @Test
    public void testDonor() {
        final MakeMetacellsFromTripletDge clp = runIt("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz",
                "donor_map.txt", Arrays.asList("DONOR"));
        TestUtils.testFilesSame(new File(TEST_DATA_DIR, "expected.donor.metacells.txt.gz"), clp.OUTPUT);
    }

    @Test
    public void testDonorCellType() {
        final MakeMetacellsFromTripletDge clp = runIt("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz",
                "donor_cell_type_map.txt", Arrays.asList("DONOR", "CELL_TYPE"));
        TestUtils.testFilesSame(new File(TEST_DATA_DIR, "expected.donor_celltype.metacells.txt.gz"), clp.OUTPUT);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testNegativeFeatureCountMismatch() {
        final MakeMetacellsFromTripletDge clp = runIt("matrix.mtx.gz", "barcodes.tsv.gz", "features_bad.tsv.gz",
                "donor_map.txt", Arrays.asList("DONOR"));
    }

    private MakeMetacellsFromTripletDge runIt(final String matrix, final String barcodes, final String features,
                                              final String mapping, final List<String> groups) {
        final MakeMetacellsFromTripletDge clp = new MakeMetacellsFromTripletDge();
        clp.TMP_DIR = Arrays.asList(TestUtils.createTempDirectory("MakeMetacellsFromTripletDgeTest."));
        clp.MATRIX = new File(TEST_DATA_DIR, matrix);
        clp.BARCODES = new File(TEST_DATA_DIR, barcodes);
        clp.FEATURES = new File(TEST_DATA_DIR, features);
        clp.MAPPING = new File(TEST_DATA_DIR, mapping);
        clp.GROUP_COLUMNS = groups;
        clp.OUTPUT = TestUtils.getTempReportFile("MakeTripletDgeTest.", ".metacells.txt.gz");
        Assert.assertEquals(clp.doWork(), 0);
        return clp;
    }
}
