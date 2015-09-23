/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2015 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
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
