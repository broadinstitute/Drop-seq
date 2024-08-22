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
package org.broadinstitute.dropseqrna.barnyard.digitalexpression;

import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.dropseqrna.barnyard.DigitalExpression;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.MergeDgeSummaries;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class MergeDgeSummariesTest {
    public static final File TESTDATA_DIR = new File("testdata/org/broadinstitute/dropseq/metrics/MergeDgeSummaries");
    public static final File GENE_SUMMARY = new File(TESTDATA_DIR, "N701.auto.exonic+intronic.digital_expression_summary.txt");
    public static final File METAGENE_SUMMARY = new File(TESTDATA_DIR, "N701.auto.exonic+intronic.metagene.digital_expression_summary.txt");
    public static final File ONE_CELL_GENE_SUMMARY = new File(TESTDATA_DIR, "one_cell.digital_expression_summary.txt");
    public static final String CBC = "CGTCCCCTCAGC";


    @Test
    public void testNonOverlapping() {
        List<DigitalExpression.DESummary> merged = doMerge(GENE_SUMMARY, ONE_CELL_GENE_SUMMARY, false);
        List<DigitalExpression.DESummary> geneSummaries = readMetrics(GENE_SUMMARY);
        DigitalExpression.DESummary oneCellSummary = readMetrics(ONE_CELL_GENE_SUMMARY).get(0);
        Assert.assertEquals(merged.size(), geneSummaries.size() + 1);
        final DigitalExpression.DESummary oneCellMerged = findMetricForCbc(merged, oneCellSummary.CELL_BARCODE);
        Assert.assertEquals(oneCellMerged, oneCellSummary);
    }

    @Test public void testOverlapping() {
        List<DigitalExpression.DESummary> mergedSummaries = doMerge(GENE_SUMMARY, METAGENE_SUMMARY, true);
        List<DigitalExpression.DESummary> geneSummaries = readMetrics(GENE_SUMMARY);
        List<DigitalExpression.DESummary> metageneSummaries = readMetrics(METAGENE_SUMMARY);
        Set<String> expectedCellBarcodes = Stream.concat(geneSummaries.stream(), metageneSummaries.stream()).map(
                deSummary -> deSummary.CELL_BARCODE).collect(Collectors.toSet());
        Assert.assertEquals(mergedSummaries.size(), expectedCellBarcodes.size());
        DigitalExpression.DESummary mergedSummary = findMetricForCbc(mergedSummaries, CBC);
        DigitalExpression.DESummary geneSummary = findMetricForCbc(geneSummaries, CBC);
        DigitalExpression.DESummary metageneSummary = findMetricForCbc(metageneSummaries, CBC);
        Assert.assertEquals(mergedSummary.NUM_GENES, geneSummary.NUM_GENES + metageneSummary.NUM_GENES);
        Assert.assertEquals(mergedSummary.NUM_GENIC_READS, geneSummary.NUM_GENIC_READS + metageneSummary.NUM_GENIC_READS);
        Assert.assertEquals(mergedSummary.NUM_TRANSCRIPTS, geneSummary.NUM_TRANSCRIPTS + metageneSummary.NUM_TRANSCRIPTS);
    }

    // There is some overlap in CBCs in the 2 inputs.  Should fail because !ACCUMULATE_CELL_BARCODE_METRICS
    @Test(expectedExceptions = RuntimeException.class)
    public void testNegative() {
        doMerge(GENE_SUMMARY, METAGENE_SUMMARY, false);
    }

    private DigitalExpression.DESummary findMetricForCbc(final List<DigitalExpression.DESummary> metrics, final String cbc) {
        final Optional<DigitalExpression.DESummary> opt = metrics.stream().filter(deSummary -> deSummary.CELL_BARCODE.equals(cbc)).findFirst();
        Assert.assertTrue(opt.isPresent());
        return opt.get();
    }

    private List<DigitalExpression.DESummary> readMetrics(final File f) {
        return MetricsFile.readBeans(f);
    }

    private List<DigitalExpression.DESummary> doMerge(final File input1, final File input2, final boolean accumulate) {
        final MergeDgeSummaries clp = new MergeDgeSummaries();
        clp.INPUT = Arrays.asList(input1, input2);
        clp.OUTPUT = TestUtils.getTempReportFile("MergeDgeSummaries.", "digital_expression.txt");
        clp.ACCUMULATE_CELL_BARCODE_METRICS = accumulate;
        Assert.assertEquals(clp.doWork(), 0);
        return readMetrics(clp.OUTPUT);
    }
}
