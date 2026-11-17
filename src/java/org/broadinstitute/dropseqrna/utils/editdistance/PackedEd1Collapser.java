/*
 * MIT License
 *
 * Copyright 2026 Broad Institute
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
package org.broadinstitute.dropseqrna.utils.editdistance;

import org.broadinstitute.dropseqrna.utils.ObjectCounter;

import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Edit distance 1 (Hamming) collapse of a group of barcodes that finds the neighbors of a barcode by looking up its
 * single-base changes in a {@link PackedBarcodeSet}, instead of comparing it against every other barcode.
 * <p>
 * Gives the same result as {@link MapBarcodesByEditDistance#collapseBarcodesGeneric} with a Hamming distance of 1:
 * the core barcodes are taken in the given order, and each absorbs the barcodes still remaining that differ from it at
 * exactly one position.  Barcodes that can't be packed (not all the same length, at most 31, or a base other than
 * upper-case ACGT) are reported by returning null from {@link #collapse}, so the caller can use the comparison search.
 * <p>
 * <b>Why there is a fast path.</b>  The comparison search takes each core barcode in turn and compares it with every
 * barcode that remains, so the work for a group of n barcodes grows with n squared.  Group sizes are very skewed: in a
 * 10x run (about 12.7 million cell/gene groups of 12 base UMIs) 94% of the groups have 50 or fewer UMIs and account
 * for about 3% of the pairwise comparisons, while the 0.3% of groups with more than 100 UMIs account for about 80%.
 * A few very large groups (highly expressed genes in large cells) therefore dominate the time of every tool that
 * collapses UMIs, such as DigitalExpression and ChimericReportEditDistanceCollapse.
 * <p>
 * For edit distance 1 there is a cheaper way to find the neighbors of a barcode: it has only 3 * length single-base
 * changes, so they can be looked up in a hash set, which takes time proportional to n instead of n squared.  The
 * barcodes are 2-bit packed into longs ({@link PackedBarcodeSet}) so that no Strings are created for the lookups.
 * <p>
 * <b>Why it is not always used.</b>  For a small group comparing a barcode with the others is a handful of
 * early-exit comparisons, which is cheaper than 3 * length hash probes, so {@link #isWorthwhile} switches to the packed
 * search only from 3 * length barcodes, where the two cost about the same number of operations.  The packed search also can't represent an N (or
 * any base other than upper-case ACGT), so such groups use the comparison search.
 * <p>
 * Collapsed barcodes are marked as removed in the packed set, so they can neither be matched again nor absorb other
 * barcodes later, as when the comparison search removes them from the list of remaining barcodes.  A randomized
 * test in {@code MapBarcodesByEditDistanceTest} checks that both searches give the same result.
 */
public class PackedEd1Collapser {
    /**
     * A barcode has 3 single-base changes per base, so finding its neighbors by lookup takes 3 * barcodeLength hash
     * probes, while comparing it against the other barcodes takes one early-exit comparison each.  The packed search
     * is used from the point where there are at least as many barcodes as lookups.
     */
    private static final int NEIGHBORS_PER_BASE = 3;

    private PackedEd1Collapser() {
    }

    /**
     * @return true if there are enough barcodes for looking up single-base changes to beat comparing every pair.  The
     * barcode length is taken from the first barcode.
     */
    public static boolean isWorthwhile(final ObjectCounter<String> barcodes) {
        if (barcodes.getSize() == 0) {
            return false;
        }
        final int barcodeLength = barcodes.getKeys().iterator().next().length();
        return barcodes.getSize() >= (long) NEIGHBORS_PER_BASE * barcodeLength;
    }

    /**
     * @param coreBarcodes the barcodes that can absorb others, in the order they do so (typically by decreasing count)
     * @param barcodes all barcodes, including the core ones
     * @return map from each core barcode that was not itself collapsed to the sorted list of barcodes collapsed into
     * it, or null if the barcodes can't be packed or a core barcode is not one of the barcodes.
     */
    public static Map<String, List<String>> collapse(final List<String> coreBarcodes, final ObjectCounter<String> barcodes) {
        final Collection<String> keys = barcodes.getKeys();
        if (keys.isEmpty()) {
            return null;
        }
        final PackedBarcodeSet remaining = PackedBarcodeSet.build(keys, keys.iterator().next().length());
        if (remaining == null) {
            return null;
        }
        final long[] corePacked = new long[coreBarcodes.size()];
        for (int i = 0; i < corePacked.length; i++) {
            corePacked[i] = PackedBarcodeSet.pack(coreBarcodes.get(i));
            if (corePacked[i] == PackedBarcodeSet.NOT_PACKABLE || remaining.get(corePacked[i]) == null) {
                return null;
            }
        }

        final Map<String, List<String>> result = new HashMap<>();
        for (int i = 0; i < corePacked.length; i++) {
            // a core that was collapsed into an earlier one is gone
            if (remaining.get(corePacked[i]) == null) {
                continue;
            }
            remaining.remove(corePacked[i]);
            final List<String> close = remaining.getEd1Matches(corePacked[i]);
            // this sorting is for consistency in output.
            Collections.sort(close);
            for (final String closeBarcode : close) {
                remaining.remove(PackedBarcodeSet.pack(closeBarcode));
            }
            result.put(coreBarcodes.get(i), close);
        }
        return result;
    }
}
