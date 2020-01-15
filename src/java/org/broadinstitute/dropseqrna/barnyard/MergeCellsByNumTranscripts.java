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
package org.broadinstitute.dropseqrna.barnyard;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.metrics.MetricsFile;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;

import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        summary = "Merges multiple CellsByNumTranscripts cell barcode lists and metrics.",
        oneLineSummary = "Merges multiple CellsByNumTranscripts cell barcode lists and metrics.",
        programGroup = DropSeq.class)
public class MergeCellsByNumTranscripts
        extends CommandLineProgram {

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "CellsByNumTranscripts reports to be merged.")
    public List<File> INPUT;

    @Argument(shortName = "IM", optional = true, doc="Optional summary metrics to merge")
    public List<File> INPUT_METRICS;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Merged CellsByNumTranscripts report")
    public File OUTPUT;

    @Argument(shortName = StandardOptionDefinitions.METRICS_FILE_SHORT_NAME, optional = true,
            doc="Optional output file for the merged summary metrics")
    public File METRICS;

    private static final Log log = Log.getInstance(MergeCellsByNumTranscripts.class);

    @Override
    protected int doWork() {
        mergeCellBarcodeFiles();
        mergeSummaryMetrics();

        return 0;
    }

    private void mergeCellBarcodeFiles() {
        Set<String> mergedBarcodes = new HashSet<>();
        for (File barcodeFile : INPUT) {
            mergedBarcodes.addAll(SelectCellsByNumTranscripts.readBarcodes(barcodeFile));
        }

        SelectCellsByNumTranscripts.writeBarcodes(OUTPUT, new ArrayList<String>(mergedBarcodes));
    }

    private void mergeSummaryMetrics() {
        if (INPUT_METRICS == null && METRICS == null) {
            return;
        } else if (INPUT_METRICS != null && METRICS == null) {
            throw new IllegalArgumentException("Metric input file list is specified, but the output metrics file is null");
        } else if (INPUT_METRICS == null) {
            throw new IllegalArgumentException("Metric input file list is null, but the output metrics file is specified");
        }

        SelectCellsByNumTranscripts.Metrics mergedMetrics = null;
        for (File reportFile : INPUT_METRICS) {
            final MetricsFile<SelectCellsByNumTranscripts.Metrics, Integer> metricsFile = getMetricsFile();
            try {
                metricsFile.read(new FileReader(reportFile));
            } catch (FileNotFoundException ex) {
                throw new RuntimeException("Error reading the metrics file " + reportFile.getAbsolutePath());
            }
            if (metricsFile.getMetrics().size() > 1) {
                throw new IllegalArgumentException("Metrics file " + reportFile.getAbsolutePath() + " contains more than one metrics");
            }
            final SelectCellsByNumTranscripts.Metrics metrics = metricsFile.getMetrics().get(0);
            if (mergedMetrics == null) {
                mergedMetrics = metrics;
            } else {
                mergedMetrics.accumulate(metrics);
            }
        }

        MetricsFile<SelectCellsByNumTranscripts.Metrics, Integer> mergedMetricsFile = getMetricsFile();
        mergedMetricsFile.addMetric(mergedMetrics);
        mergedMetricsFile.write(METRICS);
    }
}
