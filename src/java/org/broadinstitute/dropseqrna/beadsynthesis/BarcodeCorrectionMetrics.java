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
package org.broadinstitute.dropseqrna.beadsynthesis;

import htsjdk.samtools.metrics.MetricBase;

public class BarcodeCorrectionMetrics extends MetricBase {
    public long NUM_READS_EXACT_MATCH;
    public long NUM_READS_CORRECTED_SINGLE_ED1;
    public long NUM_READS_CORRECTED_MULTI_ED1;
    public long NUM_READS_UNCORRECTABLE_NO_ED1_BARCODES;
    public long NUM_READS_UNCORRECTED_AMBIGUOUS;

    public void merge(final BarcodeCorrectionMetrics that) {
        this.NUM_READS_EXACT_MATCH += that.NUM_READS_EXACT_MATCH;
        this.NUM_READS_CORRECTED_SINGLE_ED1 += that.NUM_READS_CORRECTED_SINGLE_ED1;
        this.NUM_READS_CORRECTED_MULTI_ED1 += that.NUM_READS_CORRECTED_MULTI_ED1;
        this.NUM_READS_UNCORRECTABLE_NO_ED1_BARCODES += that.NUM_READS_UNCORRECTABLE_NO_ED1_BARCODES;
        this.NUM_READS_UNCORRECTED_AMBIGUOUS += that.NUM_READS_UNCORRECTED_AMBIGUOUS;
    }
}
