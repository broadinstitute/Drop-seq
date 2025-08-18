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

import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.illumina.BarcodeMetric;
import picard.illumina.ExtractIlluminaBarcodes;

import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Merge barcode_metrics files.",
        oneLineSummary = "Merge barcode_metrics files.",
        programGroup = DropSeq.class
)
public class MergeBarcodeMetrics
        extends CommandLineProgram {
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "One or more barcode_metrics files to be merged.", minElements = 1)
    public List<File> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The file to write merged barcode_metrics to.")
    public File OUTPUT;

    public static void main(final String[] args) {
        new MergeBarcodeMetrics().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        // Load all the metrics beans
        final List<BarcodeMetric> metrics = INPUT.stream().flatMap(f -> readBeans(f).stream()).collect(Collectors.toList());

        // Remember the order of the libraries
        final Set<String> librariesInOrderSet = new LinkedHashSet<>();
        metrics.forEach(m -> librariesInOrderSet.add(m.LIBRARY_NAME));
        final List<String> librariesInOrder = new ArrayList<>(librariesInOrderSet);


        //  Group by barcode.
        final Map<String, List<BarcodeMetric>> groupedMetrics = metrics.stream().collect(Collectors.groupingBy(m -> m.BARCODE));

        // Merge each set of beans that have the same barcode
        final List<BarcodeMetric> outMetrics = groupedMetrics.values().stream().map(MergeBarcodeMetrics::combineBeans).collect(Collectors.toList());

        // Separate the no-match bean
        final BarcodeMetric noMatchMetric = outMetrics.stream().filter(m -> m.BARCODE.startsWith("N")).findFirst().get();
        outMetrics.remove(noMatchMetric);

        // Compute the values that require looking at all the beans.
        // Annoying that this method requires a Map for no good reason
        final Map<String, BarcodeMetric> outMetricMap =
                outMetrics.stream().collect(Collectors.toMap(m -> m.BARCODE, m -> m));
        ExtractIlluminaBarcodes.finalizeMetrics(outMetricMap, noMatchMetric);

        // Sort metrics based on the order the libraries appeared in the inputs.
        outMetrics.sort(Comparator.comparingInt(o -> librariesInOrder.indexOf(o.LIBRARY_NAME)));

        // Write result
        final MetricsFile<BarcodeMetric, Integer> metricsFile = getMetricsFile();
        for (final BarcodeMetric barcodeMetric : outMetrics) {
            metricsFile.addMetric(barcodeMetric);
        }
        metricsFile.addMetric(noMatchMetric);
        metricsFile.write(OUTPUT);

        return 0;
    }

    private List<BarcodeMetric> readBeans(final File file) {
        return MetricsFile.readBeans(file);
    }

    private static BarcodeMetric combineBeans(final List<BarcodeMetric> beans) {
        final BarcodeMetric outMetric = beans.get(0).copy();
        beans.forEach(outMetric::merge);
        return outMetric;
    }

}
