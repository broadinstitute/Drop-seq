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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.metrics.Header;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.TranscriptomeException;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.analysis.RnaSeqMetrics;

/**
 * @author skashin
 *
 */
@CommandLineProgramProperties(summary = "Merges rnaseq metrics files.",
        oneLineSummary = "Merges rnaseq metrics files.",
        programGroup = DropSeq.class)
public class MergeRnaSeqMetrics extends CommandLineProgram {

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input rnaseq metrics files to be merged.", minElements = 1)
    public List<File> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The merged rnaseq metrics file.")
    public File OUTPUT;


    @Override
    protected int doWork() {
        IOUtil.assertFileIsWritable(OUTPUT);

        MetricsFile<RnaSeqMetrics, Integer> mergedMetricsFile = getMetricsFile();

        Map<RnaSeqMetricsKey, List<RnaSeqMetrics>> metricsMap = new HashMap<>();
        for (File file : INPUT) {
            final MetricsFile<RnaSeqMetrics, Integer> metricsFile = new MetricsFile<>();
            metricsFile.read(IOUtil.openFileForBufferedReading(file));
            for (Header header : metricsFile.getHeaders()) {
                mergedMetricsFile.addHeader(header);
            }

            for (RnaSeqMetrics metrics : metricsFile.getMetrics()) {
                RnaSeqMetricsKey key = new RnaSeqMetricsKey(metrics.SAMPLE, metrics.LIBRARY, metrics.READ_GROUP);
                List<RnaSeqMetrics> metricsList = metricsMap.get(key);
                if (metricsList == null) {
                    metricsList = new ArrayList<RnaSeqMetrics>();
                    metricsMap.put(key, metricsList);
                }
                metricsList.add(metrics);
            }
        }

        for (List<RnaSeqMetrics> metricsList : metricsMap.values()) {
            RnaSeqMetrics metrics = (metricsList.size() == 1)
                    ? metricsList.get(0)
                    : mergeMetrics(metricsList);
            mergedMetricsFile.addMetric(metrics);
        }

        BufferedWriter writer = IOUtil.openFileForBufferedWriting(OUTPUT);
        mergedMetricsFile.write(writer);
        try {
            writer.close();
        } catch (IOException ex) {
            throw new TranscriptomeException("Problem writing file", ex);
        }

        return 0;
    }

    public RnaSeqMetrics mergeMetrics(List<RnaSeqMetrics> metricsList) {
        RnaSeqMetrics mergedMetrics = null;
        for (RnaSeqMetrics metrics : metricsList) {
            if (mergedMetrics == null) {
                mergedMetrics = metrics;
                continue;
            }

            mergedMetrics.PF_BASES += metrics.PF_BASES;
            mergedMetrics.PF_ALIGNED_BASES += metrics.PF_ALIGNED_BASES;
            if (mergedMetrics.RIBOSOMAL_BASES == null) {
                mergedMetrics.RIBOSOMAL_BASES = metrics.RIBOSOMAL_BASES;
            } else {
                mergedMetrics.RIBOSOMAL_BASES += metrics.RIBOSOMAL_BASES;
            }
            mergedMetrics.CODING_BASES += metrics.CODING_BASES;
            mergedMetrics.UTR_BASES += metrics.UTR_BASES;
            mergedMetrics.INTRONIC_BASES += metrics.INTRONIC_BASES;
            mergedMetrics.INTERGENIC_BASES += metrics.INTERGENIC_BASES;
            mergedMetrics.IGNORED_READS += metrics.IGNORED_READS;
            mergedMetrics.CORRECT_STRAND_READS += metrics.CORRECT_STRAND_READS;
            mergedMetrics.INCORRECT_STRAND_READS += metrics.INCORRECT_STRAND_READS;
            mergedMetrics.NUM_R1_TRANSCRIPT_STRAND_READS += metrics.NUM_R1_TRANSCRIPT_STRAND_READS;
            mergedMetrics.NUM_R2_TRANSCRIPT_STRAND_READS += metrics.NUM_R2_TRANSCRIPT_STRAND_READS;
            mergedMetrics.NUM_UNEXPLAINED_READS += metrics.NUM_UNEXPLAINED_READS;

            mergedMetrics.MEDIAN_CV_COVERAGE += metrics.MEDIAN_CV_COVERAGE;
            mergedMetrics.MEDIAN_5PRIME_BIAS += metrics.MEDIAN_5PRIME_BIAS;
            mergedMetrics.MEDIAN_3PRIME_BIAS += metrics.MEDIAN_3PRIME_BIAS;
            mergedMetrics.MEDIAN_5PRIME_TO_3PRIME_BIAS += metrics.MEDIAN_5PRIME_TO_3PRIME_BIAS;
        }

        double meanFactor = 1.0 / metricsList.size();
        mergedMetrics.MEDIAN_CV_COVERAGE *= meanFactor;
        mergedMetrics.MEDIAN_5PRIME_BIAS *= meanFactor;
        mergedMetrics.MEDIAN_3PRIME_BIAS *= meanFactor;
        mergedMetrics.MEDIAN_5PRIME_TO_3PRIME_BIAS *= meanFactor;

        computePercentFields(mergedMetrics);

        return(mergedMetrics);
    }

    public void computePercentFields(RnaSeqMetrics metrics) {
        if (metrics.PF_ALIGNED_BASES > 0) {
            if (metrics.RIBOSOMAL_BASES != null) {
                metrics.PCT_RIBOSOMAL_BASES =  metrics.RIBOSOMAL_BASES  / (double) metrics.PF_ALIGNED_BASES;
            }
            metrics.PCT_CODING_BASES =     metrics.CODING_BASES     / (double) metrics.PF_ALIGNED_BASES;
            metrics.PCT_UTR_BASES =        metrics.UTR_BASES        / (double) metrics.PF_ALIGNED_BASES;
            metrics.PCT_INTRONIC_BASES =   metrics.INTRONIC_BASES   / (double) metrics.PF_ALIGNED_BASES;
            metrics.PCT_INTERGENIC_BASES = metrics.INTERGENIC_BASES / (double) metrics.PF_ALIGNED_BASES;
            metrics.PCT_MRNA_BASES =       metrics.PCT_CODING_BASES + metrics.PCT_UTR_BASES;
            metrics.PCT_USABLE_BASES =     (metrics.CODING_BASES + metrics.UTR_BASES) / (double) metrics.PF_BASES;
        }

        if (metrics.CORRECT_STRAND_READS > 0 || metrics.INCORRECT_STRAND_READS > 0) {
            metrics.PCT_CORRECT_STRAND_READS = metrics.CORRECT_STRAND_READS/(double)(metrics.CORRECT_STRAND_READS + metrics.INCORRECT_STRAND_READS);
        }

        final long readsExamined = metrics.NUM_R1_TRANSCRIPT_STRAND_READS + metrics.NUM_R2_TRANSCRIPT_STRAND_READS;
        if (0 < readsExamined) {
            metrics.PCT_R1_TRANSCRIPT_STRAND_READS = metrics.NUM_R1_TRANSCRIPT_STRAND_READS / (double) readsExamined;
            metrics.PCT_R2_TRANSCRIPT_STRAND_READS = metrics.NUM_R2_TRANSCRIPT_STRAND_READS / (double) readsExamined;
        }
    }
}
