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

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import htsjdk.samtools.util.Histogram;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.TranscriptomeException;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

/**
 * @author skashin
 *
 */
@CommandLineProgramProperties(summary = "Merges read quality metrics files.",
        oneLineSummary = "Merges read quality metrics files.",
        programGroup = DropSeq.class)
public class MergeReadQualityMetrics extends CommandLineProgram {

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input read quality metrics files to be merged.", minElements = 1)
    public List<File> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The merged read quality metrics file.")
    public File OUTPUT;


    @Override
    protected int doWork() {
        IOUtil.assertFileIsWritable(OUTPUT);

        Map<String, ReadQualityMetrics> metricsMap = new TreeMap<>();
        Histogram<Integer> mergedHistogram = null;
        for (File file : INPUT) {
            MetricsFile<ReadQualityMetrics, Integer> metricsFile = new MetricsFile<>();
            metricsFile.read(IOUtil.openFileForBufferedReading(file));
            for (ReadQualityMetrics metrics : metricsFile.getMetrics()) {
                ReadQualityMetrics mergedMetrics = metricsMap.get(metrics.aggregate);
                if (mergedMetrics == null)
                    metricsMap.put(metrics.aggregate, metrics);
                else {
                    mergedMetrics.totalReads += metrics.totalReads;
                    mergedMetrics.mappedReads += metrics.mappedReads;
                    mergedMetrics.hqMappedReads += metrics.hqMappedReads;
                    mergedMetrics.hqMappedReadsNoPCRDupes += metrics.hqMappedReadsNoPCRDupes;
                }
            }

            if (mergedHistogram == null)
                mergedHistogram = metricsFile.getHistogram();
            else
                mergedHistogram.addHistogram(metricsFile.getHistogram());
        }

        MetricsFile<ReadQualityMetrics, Integer> mergedMetricsFile = new MetricsFile<>();
        // Make sure the GLOBAL metrics is added first
        mergedMetricsFile.addMetric(metricsMap.remove(GatherReadQualityMetrics.GLOBAL));
        for (ReadQualityMetrics mergedMetrics : metricsMap.values()) {
            mergedMetricsFile.addMetric(mergedMetrics);
        }

        mergedMetricsFile.setHistogram(mergedHistogram);

        BufferedWriter writer = IOUtil.openFileForBufferedWriting(OUTPUT);
        mergedMetricsFile.write(writer);
        try {
            writer.close();
        } catch (IOException ex) {
            throw new TranscriptomeException("Problem writing file", ex);
        }

        return 0;
    }
}
