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
package org.broadinstitute.dropseqrna.utils.readiterators.umisort;

import htsjdk.samtools.util.SortingCollection;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.util.ArrayList;
import java.util.List;

public class UMISortCodecTest {

    private static <T> List<T> roundTrip(final SortingCollection.Codec<T> codec, final List<T> in) {
        final ByteArrayOutputStream bos = new ByteArrayOutputStream();
        codec.setOutputStream(bos);
        in.forEach(codec::encode);
        final SortingCollection.Codec<T> reader = codec.clone();
        reader.setInputStream(new ByteArrayInputStream(bos.toByteArray()));
        final List<T> out = new ArrayList<>();
        T t;
        while ((t = reader.decode()) != null) {
            out.add(t);
        }
        return out;
    }

    private static String repeat(final String s, final int n) {
        return s.repeat(n);
    }

    @Test
    public void testUMIReadRoundTrip() {
        final List<UMIRead> in = List.of(
                new UMIRead("read1", true, "ACGT", "GENE1", "+", "CODING", "TTTT"),
                new UMIRead(null, false, "ACGT", "GENE1", null, null, "TTTT"),
                new UMIRead("", false, "", "", "", "", ""),
                new UMIRead("réad-中文", false, "AC", "GéN中", "-", "INTRONIC", "A"),
                // longer than the 64K limit of DataOutputStream.writeUTF
                new UMIRead(repeat("n", 70000), true, "C", repeat("G", 70000), null, null, "U"));
        final List<UMIRead> out = roundTrip(new UMIReadCodec(), in);
        Assert.assertEquals(out.size(), in.size());
        for (int i = 0; i < in.size(); i++) {
            Assert.assertEquals(out.get(i).toString(), in.get(i).toString());
            Assert.assertEquals(out.get(i).negativeStrand, in.get(i).negativeStrand);
            Assert.assertEquals(out.get(i).readName, in.get(i).readName);
            Assert.assertEquals(out.get(i).geneStrand, in.get(i).geneStrand);
            Assert.assertEquals(out.get(i).geneFunction, in.get(i).geneFunction);
        }
    }

    @Test
    public void testUMIKeyCountRoundTrip() {
        final List<UMIKeyCount> in = List.of(
                new UMIKeyCount("ACGT", "GENE", "TTTT", 1, 0),
                new UMIKeyCount("", "", "", 0, 5),
                new UMIKeyCount("中", "Gé", "U", 127, 127),
                new UMIKeyCount("A", "B", "C", 128, 128),
                new UMIKeyCount("A", "B", "C", Integer.MAX_VALUE, Long.MAX_VALUE),
                new UMIKeyCount(repeat("c", 70000), "g", "u", 300, 1L << 40));
        final List<UMIKeyCount> out = roundTrip(new UMIKeyCountCodec(), in);
        Assert.assertEquals(out.size(), in.size());
        for (int i = 0; i < in.size(); i++) {
            Assert.assertEquals(out.get(i).toString(), in.get(i).toString());
        }
    }

    @Test
    public void testEmptyStreamDecodesToNull() {
        Assert.assertTrue(roundTrip(new UMIReadCodec(), List.of()).isEmpty());
        Assert.assertTrue(roundTrip(new UMIKeyCountCodec(), List.of()).isEmpty());
    }

    @Test
    public void testCombineAndOrder() {
        final UMIKeyCount a = new UMIKeyCount("c", "g", "u", 2, 9);
        final UMIKeyCount b = new UMIKeyCount("c", "g", "u", 3, 4);
        Assert.assertEquals(UMIKeyCount.combine(a, b).getCount(), 5);
        Assert.assertEquals(a.getFirstSeen(), 4, "combining keeps the earlier firstSeen");
        Assert.assertThrows(ArithmeticException.class,
                () -> UMIKeyCount.combine(new UMIKeyCount("c", "g", "u", Integer.MAX_VALUE, 0), b));

        final UMIKeyCount cellA = new UMIKeyCount("a", "z", "u", 1, 0);
        final UMIKeyCount cellB = new UMIKeyCount("b", "a", "u", 1, 0);
        Assert.assertTrue(UMIKeyCount.sortOrder(true).compare(cellA, cellB) < 0);
        Assert.assertTrue(UMIKeyCount.sortOrder(false).compare(cellA, cellB) > 0);
        Assert.assertEquals(UMIKeyCount.groupOrder(true).compare(a, new UMIKeyCount("c", "g", "other", 1, 0)), 0);
    }
}
