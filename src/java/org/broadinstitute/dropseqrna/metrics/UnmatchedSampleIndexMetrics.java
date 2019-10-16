/*
 * MIT License
 *
 * Copyright 2019 Broad Institute
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

import htsjdk.samtools.metrics.MetricBase;

public class UnmatchedSampleIndexMetrics
        extends MetricBase {
    public String BARCODE;
    public int COUNT;
    public double PCT_OF_UNMATCHED_READS;
    public double PCT_OF_ALL_READS;

    public UnmatchedSampleIndexMetrics(String BARCODE, int COUNT, double PCT_OF_UNMATCHED_READS, double PCT_OF_ALL_READS) {
        this.BARCODE = BARCODE;
        this.COUNT = COUNT;
        this.PCT_OF_UNMATCHED_READS = PCT_OF_UNMATCHED_READS;
        this.PCT_OF_ALL_READS = PCT_OF_ALL_READS;
    }

    public UnmatchedSampleIndexMetrics() {
    }
}
