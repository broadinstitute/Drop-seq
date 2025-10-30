/*
 * MIT License
 *
 * Copyright 2025 Broad Institute
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

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.dropseqrna.beadsynthesis.BarcodeCorrectionMetrics;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.Reader;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class MergeBarcodeCorrectionMetricsTest {
    private static final File TEST_DATA_DIR =
            new File("testdata/org/broadinstitute/dropseq/metrics/MergeBarcodeCorrectionMetrics");
    private static final List<File> TEST_FILES =
            IntStream.rangeClosed(1, 2).mapToObj(i -> new File(TEST_DATA_DIR, String.format("%d.corrected_barcode_metrics", i))).
                    collect(Collectors.toList());
    @Test
    public void testMerge() throws Exception {
        final MergeBarcodeCorrectionMetrics mbm = new MergeBarcodeCorrectionMetrics();
        mbm.INPUT = TEST_FILES;
        final File outputMetricsFile = File.createTempFile("MergeBarcodeCorrectionMetricsTest.", ".corrected_barcode_metrics");
        outputMetricsFile.deleteOnExit();
        mbm.OUTPUT = outputMetricsFile;
        mbm.DELETE_INPUTS = false;
        Assert.assertEquals(mbm.doWork(), 0);
        final MetricsFile<BarcodeCorrectionMetrics, Integer> outputMetrics = readMetricsFile(outputMetricsFile);
        final List<MetricsFile<BarcodeCorrectionMetrics, Integer>> inputMetrics = TEST_FILES.stream().map(f -> {
            try {
                return readMetricsFile(f);
            } catch (IOException e) {
                // annoying lambda exception handling
                throw new RuntimeException(e);
            }
        }).collect(Collectors.toList());
        Assert.assertEquals(outputMetrics.getMetrics().size(), 1);
        // Spot check the metrics and histogram
        long expectedNumReadsExactMatch = inputMetrics.stream().
                mapToLong(imf -> imf.getMetrics().get(0).NUM_READS_EXACT_MATCH).
                sum();
        long expected1ED1Candidates = inputMetrics.stream().
                mapToInt(imf -> (int)imf.getHistogram().get(1).getValue()).sum();
        Assert.assertEquals(outputMetrics.getMetrics().get(0).NUM_READS_EXACT_MATCH, expectedNumReadsExactMatch);
        Assert.assertEquals((int)outputMetrics.getHistogram().get(1).getValue(), expected1ED1Candidates);
    }

    private MetricsFile<BarcodeCorrectionMetrics, Integer> readMetricsFile(File outputMetricsFile) throws IOException {
        final MetricsFile<BarcodeCorrectionMetrics, Integer> mf = new MetricsFile<>();
        final Reader reader = IOUtil.openFileForBufferedReading(outputMetricsFile);
        try {
            mf.read(reader);
        } finally {
            CloserUtil.close(reader);
        }
        return mf;
    }
}
