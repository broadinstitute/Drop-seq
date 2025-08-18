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
package org.broadinstitute.dropseqrna.metrics;

import java.io.File;
import java.util.List;
import java.util.stream.Collectors;

import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.RExecutor;

/**
 * @author skashin
 *
 */
@CommandLineProgramProperties(summary = "Merges mean quality by cycle reports.",
        oneLineSummary = "Merges mean quality by cycle reports",
        programGroup = DropSeq.class)
public class MergeMeanQualityByCycle extends CommandLineProgram {

    private static final Log log = Log.getInstance(MergeMeanQualityByCycle.class);

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input mean quality by cycle reports to be merged", minElements = 1)
    public List<File> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The merged mean quality by cycle report.")
    public File OUTPUT;

    @Argument(shortName="CHART", doc="A file (with .pdf extension) to write the chart to.")
    public File CHART_OUTPUT;


    @Override
    protected int doWork() {
        IOUtil.assertFileIsWritable(OUTPUT);

        MetricsFile<?, Integer> mergedMetricsFile = getMetricsFile();

        Histogram<Integer> mergedHistogram = null;
        for (File file : INPUT) {
            final MetricsFile<?, Integer> metricsFile = new MetricsFile<>();
            metricsFile.read(IOUtil.openFileForBufferedReading(file));
            for (Header header : metricsFile.getHeaders()) {
                mergedMetricsFile.addHeader(header);
            }

            if (mergedHistogram == null) {
                mergedHistogram = metricsFile.getHistogram();
            } else {
                mergedHistogram.addHistogram(metricsFile.getHistogram());
            }
        }

        Histogram<Integer> meanFactorHistogram = new Histogram<Integer>();
        for (Integer key : mergedHistogram.keySet()) {
            meanFactorHistogram.increment(key, INPUT.size());
        }
        mergedHistogram = mergedHistogram.divideByHistogram(meanFactorHistogram);
        mergedHistogram.setBinLabel("CYCLE");
        mergedHistogram.setValueLabel("MEAN_QUALITY");

        mergedMetricsFile.addHistogram(mergedHistogram);
        mergedMetricsFile.write(OUTPUT);

        if (CHART_OUTPUT != null)
            generateChart();

        return 0;
    }

    protected void generateChart() {
        // Get all filenames that contributed to this chart.
        final List<String> inputFilenames = INPUT.stream().map(File::getName).collect(Collectors.toList());
        // Run R to generate a chart
        final int rResult = RExecutor.executeFromClasspath(
                "picard/analysis/meanQualityByCycle.R",
                OUTPUT.getAbsolutePath(),
                CHART_OUTPUT.getAbsolutePath(),
                StringUtil.join("\n", inputFilenames),
                "");

        if (rResult != 0) {
            throw new PicardException("R script meanQualityByCycle.R failed with return code " + rResult);
        }
    }
}
