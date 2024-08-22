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
package org.broadinstitute.dropseqrna.metagene;

import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;


public class MergeMetaGeneReportsTest {
    static final File SMN_FILE = new File("testdata/org/broadinstitute/dropseq/metagene/BA46_downsampled_SMN.bam");

    @Test
    public void testDoWork() {
        File expectedMetaGenesReport = TestUtils.getTempReportFile("MergeMetaGeneReports", ".txt.gz");
        File metaGenesReport1 = TestUtils.getTempReportFile("MergeMetaGeneReports", ".txt.gz");
        File metaGenesReport2 = TestUtils.getTempReportFile("MergeMetaGeneReports", ".txt.gz");
        File mergedMetaGenesReport = TestUtils.getTempReportFile("MergeMetaGeneReports", ".txt.gz");

        List<File> splitBAMFileList = TestUtils.splitBamFile(SMN_FILE, 2);
        Assert.assertEquals(splitBAMFileList.size(), 2);

        discoverMetaGenes(SMN_FILE, expectedMetaGenesReport);
        discoverMetaGenes(splitBAMFileList.get(0), metaGenesReport1);
        discoverMetaGenes(splitBAMFileList.get(1), metaGenesReport2);

        MergeMetaGeneReports merger = new MergeMetaGeneReports();
        merger.INPUT = Arrays.asList(metaGenesReport1, metaGenesReport2);
        merger.OUTPUT = mergedMetaGenesReport;
        Assert.assertEquals(merger.doWork(), 0);

        UMIMetaGeneAggregation expectedAggregation = DiscoverMetaGenes.readReport(expectedMetaGenesReport);
        UMIMetaGeneAggregation mergedAggregation = DiscoverMetaGenes.readReport(mergedMetaGenesReport);

        Assert.assertTrue(mergedAggregation.equals(expectedAggregation));
    }

    private void discoverMetaGenes(final File inputBAM, final File metaGenesReport) {
        DiscoverMetaGenes discoverer = new DiscoverMetaGenes();
        discoverer.INPUT = inputBAM;
        discoverer.MIN_READ_MQ = 2;
        discoverer.META_GENE_RATIO = 0d;
        discoverer.REPORT = metaGenesReport;

        Assert.assertEquals(discoverer.doWork(), 0);
    }
}
