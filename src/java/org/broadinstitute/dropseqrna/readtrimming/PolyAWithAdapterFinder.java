/*
 * MIT License
 *
 * Copyright 2017 Broad Institute
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
package org.broadinstitute.dropseqrna.readtrimming;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.StringUtil;
import picard.util.ClippingUtility;

public class PolyAWithAdapterFinder implements PolyAFinder {

    private static char A = 'A';

    private final AdapterDescriptor adapter;
    private final int minAdapterMatch;
    private final double maxAdapterErrorRate;
    private final int minPolyALength;
    private final int minPolyALengthNoAdapterMatch;
    private final double maxPolyAErrorRate;
    private final int dubiousAdapterMatchLength;

    public PolyAWithAdapterFinder(
            AdapterDescriptor adapter,
            int minAdapterMatch,
            double maxAdapterErrorRate,
            int minPolyALength,
            int minPolyALengthNoAdapterMatch,
            double maxPolyAErrorRate,
            int dubiousAdapterMatchLength) {
        this.adapter = adapter;
        this.minAdapterMatch = minAdapterMatch;
        this.maxAdapterErrorRate = maxAdapterErrorRate;
        this.minPolyALength = minPolyALength;
        this.minPolyALengthNoAdapterMatch = minPolyALengthNoAdapterMatch;
        this.maxPolyAErrorRate = maxPolyAErrorRate;
        this.dubiousAdapterMatchLength = dubiousAdapterMatchLength;
    }

    @Override
    public PolyARun getPolyAStart(final SAMRecord r) {
        return getPolyAStart(r.getReadString(), adapter.getAdapterSequence(r));
    }

    public PolyARun getPolyAStart(final String readString, final String adapterSequence) {
        final byte[] readBases = StringUtil.stringToBytes(readString);
        int adapterClipPosition = ClippingUtility.findIndexOfClipSequence(
                readBases,
                StringUtil.stringToBytes(adapterSequence),
                minAdapterMatch,
                maxAdapterErrorRate);
        if (adapterClipPosition == ClippingUtility.NO_MATCH) {
            adapterClipPosition = readString.length();
        } else if (adapterClipPosition == 0) {
            return new PolyARun(0, 0, 0);
        }
        final SimplePolyAFinder.PolyARun ret = getPolyARun(readString, adapterClipPosition);

        // If there was a short adapter match, but not enough poly A before it,
        // see if there would be enough poly A if the adapter considered not to match.
        if (ret.isNoMatch() && adapterClipPosition < readString.length() &&
                adapterClipPosition + dubiousAdapterMatchLength >= readString.length()) {
            // If did not find enough polyA looking before adapter, try again looking from end of read.
            final SimplePolyAFinder.PolyARun tryWithoutAdapter = getPolyARun(readString, readString.length());
            if (!tryWithoutAdapter.isNoMatch()) {
                return tryWithoutAdapter;
            }
        }
        return ret;

    }

    private PolyAFinder.PolyARun getPolyARun(String readBases, int adapterClipPosition) {
        // Note whether there was actually adapter found, as opposed to just starting at the end of the read.
        final int realAdapterClipPosition;
        if (adapterClipPosition < readBases.length() - 1) {
            realAdapterClipPosition = adapterClipPosition;
        } else {
            realAdapterClipPosition = -1;
        }

        final int minPolyABases;
        if (realAdapterClipPosition == -1) {
            // If no adapter match, then allow a short poly A match because most of the poly A tail could be
            // beyond the end of the read.
            minPolyABases = minPolyALengthNoAdapterMatch;
        } else {
            // If adapter match, then there should be at least MIN_POLY_A_LENGTH poly A match, or
            // the number of bases before the adapter starts, whichever is shorter.
            minPolyABases = Math.min(adapterClipPosition, minPolyALength);
        }
        int numMisMatches = 0;

        int bestPolyARunStart = -1;
        double bestErrorRate = 1.0;

        // Start just before the adapter
        for (int i = adapterClipPosition - 1; i >= 0; --i) {
            if (readBases.charAt(i) == A) {
                final double errorRate = numMisMatches/(double)(adapterClipPosition - i);
                if (adapterClipPosition - i >= minPolyABases && errorRate <= maxPolyAErrorRate) {
                    bestPolyARunStart = i;
                    bestErrorRate = errorRate;
                }
            } else {
                ++numMisMatches;
            }
        }
        if (bestErrorRate <= maxPolyAErrorRate && adapterClipPosition - bestPolyARunStart >= minPolyABases) {
            return new SimplePolyAFinder.PolyARun(bestPolyARunStart, adapterClipPosition - bestPolyARunStart, realAdapterClipPosition);
        } else {
            return new SimplePolyAFinder.PolyARun(SimplePolyAFinder.NO_MATCH, 0, realAdapterClipPosition);
        }

    }

}
