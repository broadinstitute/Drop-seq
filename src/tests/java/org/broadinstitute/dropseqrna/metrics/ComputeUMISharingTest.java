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

import org.broadinstitute.dropseqrna.barnyard.GeneFunctionCommandLineBase;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

public class ComputeUMISharingTest {
    public static final String MAPPED = "mapped";
    public static final String UNMAPPED = "unmapped";
    private static File TESTDATA_DIR = new File("testdata/org/broadinstitute/dropseq/metrics");
    private static String INPUT_BAM = "compute_umi_sharing.%s.sam";
    private static String EXPECTED_SINGLE_COUNT_TAG_METRICS = "compute_umi_sharing.single_count_tag.%s.%d.umi_sharing_metrics";
    private static String EXPECTED_MULTI_COUNT_TAG_METRICS = "compute_umi_sharing.multi_count_tag.%d.%s.umi_sharing_metrics";

    @Test(dataProvider = "singleCountTagDataProvider")
    public void testSingleCountTag(final int editDistance, final boolean mapped) throws IOException {
        final String mappedStr;
        if (mapped) {
            mappedStr = MAPPED;
        } else {
            mappedStr = UNMAPPED;
        }
        final File outFile = File.createTempFile("ComputeUMISharing.ED" + editDistance + "." + mappedStr + ".", ".edit_distance_metrics");
        outFile.deleteOnExit();
        final ComputeUMISharing clp = new ComputeUMISharing();
        clp.EDIT_DISTANCE = Collections.singletonList(editDistance);
        clp.COUNT_TAG = Collections.singletonList("XM");
        clp.OUTPUT = outFile;
        clp.COLLAPSE_TAG = (mapped? "ZC": "rm");
        clp.UNCOLLAPSED_TAG = (mapped? "XC": "rb");
        clp.FIND_INDELS = false;
        clp.INPUT =  new File(TESTDATA_DIR, String.format(INPUT_BAM, mappedStr));
        clp.NUM_THREADS = 2;
        if (mapped) {
            clp.LOCUS_FUNCTION_LIST = GeneFunctionCommandLineBase.DEFAULT_LOCUS_FUNCTION_LIST;
        }
        Assert.assertEquals(clp.doWork(), 0);
        final File expectedMetricsFile = new File(TESTDATA_DIR, String.format(EXPECTED_SINGLE_COUNT_TAG_METRICS, mappedStr, editDistance));
        Assert.assertTrue(TestUtils.testMetricsFilesEqual(expectedMetricsFile, outFile),
                String.format("%s and %s differ", expectedMetricsFile.getAbsolutePath(), outFile.getAbsolutePath()));
    }

    @DataProvider(name="singleCountTagDataProvider")
    public Object[][] singleCountTagDataProvider() {
        ArrayList<Object[]> ret = new ArrayList<>();
        int[] editDistances = {0,1};
        boolean[] mappeds = {false, true};
        for (int editDistance : editDistances) {
            for (boolean mapped : mappeds) {
                ret.add(new Object[]{editDistance, mapped});
            }
        }
        return ret.toArray(new Object[ret.size()][]);
    }

    @Test(dataProvider = "multipleCountTagDataProvider")
    public void testMultipleCountTag(final int editDistance, final boolean geneMatch) throws IOException {
        final File outFile = File.createTempFile("ComputeUMISharingMultiTag.ED" + editDistance + "." , ".edit_distance_metrics");
        outFile.deleteOnExit();
        final ComputeUMISharing clp = new ComputeUMISharing();
        clp.EDIT_DISTANCE = Collections.singletonList(editDistance);
        clp.COUNT_TAG = (geneMatch? Arrays.asList("XM", GeneFunctionCommandLineBase.DEFAULT_GENE_NAME_TAG, GeneFunctionCommandLineBase.DEFAULT_GENE_STRAND_TAG):
                Arrays.asList("XM", GeneFunctionCommandLineBase.DEFAULT_GENE_STRAND_TAG));
        clp.OUTPUT = outFile;
        clp.COLLAPSE_TAG = "ZC";
        clp.UNCOLLAPSED_TAG = "XC";
        clp.FIND_INDELS = false;
        clp.INPUT =  new File(TESTDATA_DIR, String.format(INPUT_BAM, MAPPED));
        clp.NUM_THREADS = 2;
        clp.LOCUS_FUNCTION_LIST = GeneFunctionCommandLineBase.DEFAULT_LOCUS_FUNCTION_LIST;
        Assert.assertEquals(clp.doWork(), 0);
        final File expectedMetricsFile = new File(TESTDATA_DIR, String.format(EXPECTED_MULTI_COUNT_TAG_METRICS, editDistance, geneMatch));
        Assert.assertTrue(TestUtils.testMetricsFilesEqual(expectedMetricsFile, outFile),
                String.format("%s and %s differ", expectedMetricsFile.getAbsolutePath(), outFile.getAbsolutePath()));

    }

    @DataProvider(name="multipleCountTagDataProvider")
    public Object[][] multipleCountTagDataProvider() {
        ArrayList<Object[]> ret = new ArrayList<>();
        int[] editDistances = {0,1};
        boolean[] geneMatches = {false, true};
        for (int editDistance : editDistances) {
            for (boolean geneMatch : geneMatches) {
                ret.add(new Object[]{editDistance, geneMatch});
            }
        }
        return ret.toArray(new Object[ret.size()][]);
    }
}
