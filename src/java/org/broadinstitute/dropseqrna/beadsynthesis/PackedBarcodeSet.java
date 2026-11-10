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
package org.broadinstitute.dropseqrna.beadsynthesis;

import org.broadinstitute.dropseqrna.utils.BaseRange;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

/**
 * A set of equal-length ACGT-only barcodes, stored as 2-bit-packed longs in an open-addressing table, so that a barcode
 * can be looked up straight from read bases, and its edit-distance-1 neighbors can be looked up, without creating any
 * Strings.  Each lookup returns the String that was passed to {@link #build}, which is typically interned by the caller.
 * <p>
 * Used by {@link BarcodeCorrector} as a fast path.  A barcode that cannot be packed (a base other than upper-case ACGT,
 * or a range that runs past the end of the read) is reported as {@link #NOT_PACKABLE}, and the caller falls back to
 * its String-based logic, so the results are the same as that logic gives.
 */
class PackedBarcodeSet {
    /** Returned by {@link #packRead} when the barcode cannot be packed. */
    static final long NOT_PACKABLE = -1;
    /** 2 bits per base in a long, leaving the top bit free so packed values are never negative. */
    static final int MAX_LENGTH = 31;

    private static final byte[] BASE_CODE = new byte[256];
    static {
        Arrays.fill(BASE_CODE, (byte) -1);
        BASE_CODE['A'] = 0;
        BASE_CODE['C'] = 1;
        BASE_CODE['G'] = 2;
        BASE_CODE['T'] = 3;
    }

    /**
     * For each base code, the codes of the three other bases, in the same order as the {@code gct, agt, act, acg}
     * arrays in BarcodeCorrector, so that this search visits candidates in the same order as the String-based one.
     * That order almost never changes a result (a tie for the best likelihood cannot reach the default
     * LIKELIHOOD_RATIO, and summation order only matters at the last bit), but keeping it makes the two searches
     * identical by construction.
     */
    private static final int[][] OTHER_BASE_CODES = {
            {2, 1, 3}, // A -> G, C, T
            {0, 2, 3}, // C -> A, G, T
            {0, 1, 3}, // G -> A, C, T
            {0, 1, 2}  // T -> A, C, G
    };

    private final int barcodeLength;
    // packed barcode + 1, so that 0 can mean "empty slot"
    private final long[] keys;
    private final String[] values;
    private final int mask;
    private final int hashShift;

    private PackedBarcodeSet(final int barcodeLength, final int capacityBits) {
        this.barcodeLength = barcodeLength;
        final int capacity = 1 << capacityBits;
        this.keys = new long[capacity];
        this.values = new String[capacity];
        this.mask = capacity - 1;
        this.hashShift = 64 - capacityBits;
    }

    /**
     * @param barcodes distinct barcodes
     * @param barcodeLength the length that every barcode must have
     * @return a set, or null if the fast path cannot represent these barcodes: length out of range, a barcode of a
     * different length, or a base other than upper-case ACGT.
     */
    static PackedBarcodeSet build(final Collection<String> barcodes, final int barcodeLength) {
        if (barcodeLength < 1 || barcodeLength > MAX_LENGTH) {
            return null;
        }
        int capacityBits = 4;
        // keep the load factor at or below 0.5
        while ((1L << capacityBits) < 2L * barcodes.size()) {
            ++capacityBits;
        }
        final PackedBarcodeSet ret = new PackedBarcodeSet(barcodeLength, capacityBits);
        for (final String barcode : barcodes) {
            if (barcode.length() != barcodeLength) {
                return null;
            }
            long packed = 0;
            for (int i = 0; i < barcodeLength; ++i) {
                final char c = barcode.charAt(i);
                final int code = c < BASE_CODE.length ? BASE_CODE[c] : -1;
                if (code < 0) {
                    return null;
                }
                packed = (packed << 2) | code;
            }
            ret.put(packed, barcode);
        }
        return ret;
    }

    /**
     * Pack the barcode that is the concatenation of the ranges of the bases.
     * @param bases read bases, as in {@link htsjdk.samtools.SAMRecord#getReadBases()}
     * @return the packed barcode, or {@link #NOT_PACKABLE} if a range runs past the end of the bases or a base is not
     * upper-case ACGT.
     */
    static long packRead(final List<BaseRange> ranges, final byte[] bases) {
        long packed = 0;
        for (final BaseRange range : ranges) {
            if (range.getEnd() > bases.length) {
                return NOT_PACKABLE;
            }
            for (int i = range.getStart() - 1; i < range.getEnd(); ++i) {
                final int code = BASE_CODE[bases[i] & 0xff];
                if (code < 0) {
                    return NOT_PACKABLE;
                }
                packed = (packed << 2) | code;
            }
        }
        return packed;
    }

    /** @return the member that was packed as this value, or null if there is none. */
    String get(final long packed) {
        final long key = packed + 1;
        int idx = slot(key);
        while (true) {
            final long k = keys[idx];
            if (k == 0) {
                return null;
            }
            if (k == key) {
                return values[idx];
            }
            idx = (idx + 1) & mask;
        }
    }

    /**
     * @return the members that differ from this barcode at exactly one position, ordered by position and then by the
     * order of the other three bases, as BarcodeCorrector's String-based search orders them.
     */
    List<String> getEd1Matches(final long packed) {
        final List<String> ret = new ArrayList<>();
        for (int i = 0; i < barcodeLength; ++i) {
            final int shift = 2 * (barcodeLength - 1 - i);
            final int original = (int) ((packed >>> shift) & 3);
            final long cleared = packed & ~(3L << shift);
            for (final int other : OTHER_BASE_CODES[original]) {
                final String match = get(cleared | ((long) other << shift));
                if (match != null) {
                    ret.add(match);
                }
            }
        }
        return ret;
    }

    private void put(final long packed, final String barcode) {
        final long key = packed + 1;
        int idx = slot(key);
        while (keys[idx] != 0) {
            idx = (idx + 1) & mask;
        }
        keys[idx] = key;
        values[idx] = barcode;
    }

    private int slot(final long key) {
        return (int) ((key * 0x9E3779B97F4A7C15L) >>> hashShift);
    }
}
