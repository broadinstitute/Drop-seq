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
package org.broadinstitute.dropseqrna.metrics;

import com.google.common.collect.ImmutableMap;
import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.illumina.BarcodeMetric;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class MergeBarcodeMetricsTest {

    private static final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/metrics/MergeBarcodeMetrics");
    private static final List<File> TEST_FILES =
            IntStream.rangeClosed(1, 4).mapToObj(i -> new File(TEST_DATA_DIR, String.format("HLMCKBGXY.%d.barcode_metrics", i))).
                    collect(Collectors.toList());

    private static final Map<String, Long> PF_READS = ImmutableMap.<String, Long>builder()
            .put("AGGCAGAA", 1238826L)
            .put("CGTACTAG", 1101456L)
            .put("GGACTCCT", 2492003L)
            .put("NNNNNNNN", 391878666L)
            .put("TAAGGCGA", 1265209L)
            .put("TCCTGAGC", 1030185L).build();
    private static final long TOTAL_PF_READS = PF_READS.values().stream().mapToLong(l -> l).sum();

    @Test
    public void testMerge() throws Exception{
        final MergeBarcodeMetrics mbm = new MergeBarcodeMetrics();
        mbm.INPUT = TEST_FILES;
        final File outputMetricsFile = File.createTempFile("MergeBarcodeMetricsTest.", ".barcode_metrics");
        outputMetricsFile.deleteOnExit();
        mbm.OUTPUT = outputMetricsFile;
        Assert.assertEquals(mbm.doWork(), 0);
        final List<BarcodeMetric> beans = MetricsFile.readBeans(outputMetricsFile);
        Assert.assertEquals(beans.size(), PF_READS.size());
        for (final BarcodeMetric bean: beans) {
            final long expectedPfReads = PF_READS.get(bean.BARCODE);
            Assert.assertEquals(bean.PF_READS, expectedPfReads);
            final double expectedPct = expectedPfReads / (double) TOTAL_PF_READS;
            Assert.assertEquals(bean.PF_PCT_MATCHES, expectedPct, expectedPct/1000);
        }
    }
}
