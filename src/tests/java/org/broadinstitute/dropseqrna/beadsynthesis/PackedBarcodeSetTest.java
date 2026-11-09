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
import org.testng.Assert;
import org.testng.annotations.Test;

import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

public class PackedBarcodeSetTest {
    private static final char[] BASES = {'A', 'C', 'G', 'T'};
    // Same order as the arrays in BarcodeCorrector: the other bases for A, C, G, T
    private static final String[] OTHERS = {"GCT", "AGT", "ACT", "ACG"};

    private static String randomBarcode(final Random random, final int length) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < length; ++i) {
            sb.append(BASES[random.nextInt(4)]);
        }
        return sb.toString();
    }

    private static byte[] bytes(final String s) {
        return s.getBytes(StandardCharsets.ISO_8859_1);
    }

    /** Straightforward String implementation of the ED1 search, in the order the corrector has always used. */
    private static List<String> referenceEd1Matches(final String barcode, final Set<String> allowed) {
        final List<String> ret = new ArrayList<>();
        for (int i = 0; i < barcode.length(); ++i) {
            final int original = "ACGT".indexOf(barcode.charAt(i));
            for (final char other : OTHERS[original].toCharArray()) {
                final String candidate = barcode.substring(0, i) + other + barcode.substring(i + 1);
                if (allowed.contains(candidate)) {
                    ret.add(candidate);
                }
            }
        }
        return ret;
    }

    @Test
    public void testLookupAndEd1MatchesAgreeWithReference() {
        final Random random = new Random(1);
        // Short barcodes, so that there are many ED1 neighbors among the members.
        for (final int length : new int[]{1, 2, 8, 16, PackedBarcodeSet.MAX_LENGTH}) {
            final Set<String> members = new LinkedHashSet<>();
            final int numMembers = Math.min(2000, (int) Math.min(Integer.MAX_VALUE, 3L << (2 * Math.min(length, 10))) / 4);
            while (members.size() < numMembers) {
                members.add(randomBarcode(random, length));
            }
            final PackedBarcodeSet set = PackedBarcodeSet.build(members, length);
            Assert.assertNotNull(set);
            final List<BaseRange> ranges = BaseRange.parseBaseRange("1-" + length);
            for (int i = 0; i < 3000; ++i) {
                final String query = random.nextBoolean() ?
                        members.stream().skip(random.nextInt(members.size())).findFirst().get() :
                        randomBarcode(random, length);
                final long packed = PackedBarcodeSet.packRead(ranges, bytes(query));
                Assert.assertNotEquals(packed, PackedBarcodeSet.NOT_PACKABLE);
                Assert.assertEquals(set.get(packed), members.contains(query) ? query : null);
                Assert.assertEquals(set.getEd1Matches(packed), referenceEd1Matches(query, members), query);
            }
        }
    }

    @Test
    public void testAllZerosBarcodeIsNotConfusedWithEmptySlot() {
        final PackedBarcodeSet set = PackedBarcodeSet.build(List.of("AAAA", "AAAC"), 4);
        Assert.assertNotNull(set);
        final List<BaseRange> ranges = BaseRange.parseBaseRange("1-4");
        Assert.assertEquals(set.get(PackedBarcodeSet.packRead(ranges, bytes("AAAA"))), "AAAA");
        Assert.assertNull(set.get(PackedBarcodeSet.packRead(ranges, bytes("AAAG"))));
        Assert.assertEquals(set.getEd1Matches(PackedBarcodeSet.packRead(ranges, bytes("AAAG"))), List.of("AAAC", "AAAA").reversed());
    }

    @Test
    public void testBuildReturnsNullWhenBarcodesCannotBePacked() {
        Assert.assertNull(PackedBarcodeSet.build(List.of("ACGT", "ACGN"), 4));
        Assert.assertNull(PackedBarcodeSet.build(List.of("ACGT", "acgt"), 4));
        Assert.assertNull(PackedBarcodeSet.build(List.of("ACGT", "ACG"), 4));
        Assert.assertNull(PackedBarcodeSet.build(List.of("ACGT"), PackedBarcodeSet.MAX_LENGTH + 1));
        Assert.assertNull(PackedBarcodeSet.build(List.of(""), 0));
        Assert.assertNotNull(PackedBarcodeSet.build(List.of(), 4));
    }

    @Test
    public void testPackReadNotPackable() {
        final List<BaseRange> ranges = BaseRange.parseBaseRange("1-4");
        Assert.assertEquals(PackedBarcodeSet.packRead(ranges, bytes("ACGN")), PackedBarcodeSet.NOT_PACKABLE);
        Assert.assertEquals(PackedBarcodeSet.packRead(ranges, bytes("acgt")), PackedBarcodeSet.NOT_PACKABLE);
        Assert.assertEquals(PackedBarcodeSet.packRead(ranges, bytes("ACG")), PackedBarcodeSet.NOT_PACKABLE);
        // bases outside of the ranges do not matter
        Assert.assertNotEquals(PackedBarcodeSet.packRead(ranges, bytes("ACGTNNNN")), PackedBarcodeSet.NOT_PACKABLE);
    }

    @Test
    public void testMultipleRanges() {
        final List<BaseRange> ranges = BaseRange.parseBaseRange("1-2:5-6");
        final PackedBarcodeSet set = PackedBarcodeSet.build(List.of("ACGT"), 4);
        Assert.assertNotNull(set);
        // positions 1,2,5,6 of this read are A,C,G,T
        Assert.assertEquals(set.get(PackedBarcodeSet.packRead(ranges, bytes("ACTTGTNN"))), "ACGT");
        Assert.assertNull(set.get(PackedBarcodeSet.packRead(ranges, bytes("ACGTTTNN"))));
    }
}
