/*
 * MIT License
 *
 * Copyright 2020 Broad Institute
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

import org.broadinstitute.dropseqrna.eqtl.PrepareEqtlSnpGeneMap;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.File;
import java.util.HashSet;
import java.util.Set;

public class PrepareEqtlSnpGeneMapTest {
    private final File TEST_DATA_DIR =
            new File("testdata/org/broadinstitute/dropseq/eqtl/PrepareEqtlSnpGeneMap");

    private final File NON_SPECIFIC_INTERVALS =
            new File(TEST_DATA_DIR, "ATAC-Seq.GSM2264802_C15.GRCh38.chrX.subset.interval_list");
    private final File GENE_SPECIFIC_INTERVALS =
            new File(TEST_DATA_DIR, "Hi-C.iPSC-1.GRCh38.chrX.subset.interval_list");
    private final File SNP_LOCATIONS = new File(TEST_DATA_DIR, "chrX.subset.variant_locations.txt");
    private final File GENE_LOCATIONS = new File(TEST_DATA_DIR, "chrX.gene_locations.txt");

    private final int CIS_DIST = 9824;
    private final int NS_CIS_DIST = 293659;
    private final int GS_CIS_DIST = 260123;

    // These are not the only results, but they will drop out
    // when the given distances are reduced.
    private final String[] EXPECTED_MAPPINGS = {
            "chrX:266498:C:A PLCXD1", // cisDist 9824
            "chrX:54547476:G:A MAGED2", // gsCisDist 260123
            "chrX:54639302:T:C PFKFB1" // nsCisDist 293659
    };

    private Set<String> prepareHelper(final int cisDist, final int nsCisDist, final int gsCisDist) {
        final PrepareEqtlSnpGeneMap clp = new PrepareEqtlSnpGeneMap();
        clp.SNP_LOCATIONS = SNP_LOCATIONS;
        clp.GENE_LOCATION_FILE = GENE_LOCATIONS;
        clp.CIS_DIST = cisDist;
        clp.NS_INTERVAL_LIST = NON_SPECIFIC_INTERVALS;
        clp.NS_CIS_DIST = nsCisDist;
        clp.GS_INTERVAL_LIST = GENE_SPECIFIC_INTERVALS;
        clp.GS_CIS_DIST = gsCisDist;
        clp.OUTPUT = TestUtils.getTempReportFile("snp_gene_map.", ".txt");
        Assert.assertEquals(clp.doWork(), 0);
        final TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(clp.OUTPUT);
        final Set<String> ret = new HashSet<>();
        for (final TabbedTextFileWithHeaderParser.Row row: parser) {
            ret.add(
               row.getField("snp") + " " + row.getField("gene")
            );
        }
        return ret;
    }

    @Test
    public void testPositive() {
        final Set<String> mappings = prepareHelper(CIS_DIST, NS_CIS_DIST, GS_CIS_DIST);
        for (final String expected: EXPECTED_MAPPINGS) {
            Assert.assertTrue(mappings.contains(expected));
        }
    }

    // Confirm that when the distances are reduced by 1, the mapping disappear.
    @Test
    public void testNegative() {
        final Set<String> mappings = prepareHelper(CIS_DIST-1, NS_CIS_DIST-1, GS_CIS_DIST-1);
        for (final String expected: EXPECTED_MAPPINGS) {
            Assert.assertFalse(mappings.contains(expected));
        }
    }


}
