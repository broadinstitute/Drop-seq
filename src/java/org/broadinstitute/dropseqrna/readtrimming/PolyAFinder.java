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

public interface PolyAFinder {
    public static int NO_MATCH = -1;

    /**
	 * Gets the first base of a polyA run in the sequence, or -1 if there is no polyA run of at least minNumBases.
	 * @param rec the read to test.
	 * @return The index of the first base of the polyA run.  This is a 0 based index!
	 */
    PolyARun getPolyAStart(final SAMRecord rec);

    public static class PolyARun {
        /** 0-based start of poly A run, or NO_MATCH if none found */
        public final int startPos;
        public final int length;
        /** 0-based start of adapter, or NO_MATCH if none found */
        public final int adapterStartPos;

        public PolyARun(int startPos, int length) {
            this(startPos, length, -1);
        }

        public PolyARun(int startPos, int length, int adapterStartPos) {
            this.startPos = startPos;
            this.length = length;
            this.adapterStartPos = adapterStartPos;
        }

        // Inclusive end position
        public int endPos() {return startPos + length - 1; }

        boolean isNoMatch() {
            return this.startPos == NO_MATCH;
        }
        public static PolyARun NO_MATCH_RUN = new PolyARun(NO_MATCH, 0);

    }
}
