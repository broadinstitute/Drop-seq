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
package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;

import java.util.function.Predicate;

/**
 * true if SAMRecord passes filter
 */
public class MapQualityPredicate
        implements Predicate<SAMRecord> {

    private final Integer mapQuality;
    private final boolean rejectNonPrimaryReads;

    public MapQualityPredicate(Integer mapQuality, boolean rejectNonPrimaryReads) {
        this.mapQuality = mapQuality;
        this.rejectNonPrimaryReads = rejectNonPrimaryReads;
    }

    @Override
    public boolean test(SAMRecord r) {
        boolean flag=(rejectNonPrimaryReads && r.isSecondaryOrSupplementary()) ||
                (this.mapQuality!=null && r.getMappingQuality() < this.mapQuality);
        return !flag;
    }
}
