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

import java.util.List;

import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.RnaSeqMtMetrics;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;

import htsjdk.samtools.util.Log;
import picard.analysis.RnaSeqMetrics;

/**
 * @author skashin
 *
 */
@CommandLineProgramProperties(summary = "Merges single cell rnaseq metrics files",
        oneLineSummary = "Merges single cell rnaseq metrics files",
        programGroup = DropSeq.class)
public class MergeSingleCellRnaSeqMetrics extends MergeRnaSeqMetrics {

    private static final Log log = Log.getInstance(MergeSingleCellRnaSeqMetrics.class);

    @Override
    public RnaSeqMetrics mergeMetrics(List<RnaSeqMetrics> metricsList) {
        RnaSeqMetrics mergedMetrics = super.mergeMetrics(metricsList);
        if (!(mergedMetrics instanceof RnaSeqMtMetrics)) {
            throw new RuntimeException("The merged metrics is not of class RnaSeqMtMetrics");
        }

        RnaSeqMtMetrics mergedMtMetrics = (RnaSeqMtMetrics) mergedMetrics;
        for (RnaSeqMetrics metrics : metricsList) {
            mergedMtMetrics.MT_BASES += ((RnaSeqMtMetrics) metrics).MT_BASES;
        }

        if (mergedMtMetrics.PCT_RIBOSOMAL_BASES == null)
            mergedMtMetrics.PCT_RIBOSOMAL_BASES = 0d;
        if (mergedMtMetrics.PF_ALIGNED_BASES > 0)
            mergedMtMetrics.PCT_MT_BASES = mergedMtMetrics.MT_BASES / (double) mergedMtMetrics.PF_ALIGNED_BASES;

        return mergedMtMetrics;
    }
}
