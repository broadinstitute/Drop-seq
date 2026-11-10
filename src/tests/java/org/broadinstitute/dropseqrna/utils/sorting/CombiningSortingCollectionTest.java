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
package org.broadinstitute.dropseqrna.utils.sorting;

import htsjdk.samtools.util.CloseableIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.umisort.UMIKeyCount;
import org.broadinstitute.dropseqrna.utils.readiterators.umisort.UMIKeyCountCodec;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.TreeMap;

public class CombiningSortingCollectionTest {

    private static List<UMIKeyCount> run(final List<UMIKeyCount> input, final int maxInRam, final int expectedSpills) {
        final CombiningSortingCollection<UMIKeyCount> c = new CombiningSortingCollection<>(new UMIKeyCountCodec(),
                UMIKeyCount.sortOrder(true), UMIKeyCount::combine, maxInRam);
        input.forEach(c::add);
        Assert.assertEquals(c.getNumSpillFiles(), expectedSpills);
        final List<UMIKeyCount> out = new ArrayList<>();
        try (CloseableIterator<UMIKeyCount> it = c.iterator()) {
            it.forEachRemaining(out::add);
        }
        return out;
    }

    private static TreeMap<String, Integer> expected(final List<UMIKeyCount> input) {
        final TreeMap<String, Integer> m = new TreeMap<>();
        for (final UMIKeyCount k : input) {
            m.merge(k.cell + "|" + k.gene + "|" + k.umi, k.getCount(), Integer::sum);
        }
        return m;
    }

    private static List<UMIKeyCount> randomInput(final int n, final int distinct, final long seed) {
        final Random r = new Random(seed);
        final List<UMIKeyCount> in = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            final int k = r.nextInt(distinct);
            in.add(new UMIKeyCount("c" + k % 7, "g" + k % 5, "u" + k, 1 + r.nextInt(3), i));
        }
        return in;
    }

    /** The combiner merges counts in place, so the expected counts must be computed before the input is run. */
    private static void check(final List<UMIKeyCount> input, final int maxInRam, final int expectedSpills) {
        final TreeMap<String, Integer> exp = expected(input);
        final List<UMIKeyCount> out = run(input, maxInRam, expectedSpills);
        Assert.assertEquals(out.size(), exp.size(), "records must be unique");
        final List<String> expKeys = new ArrayList<>(exp.keySet());
        // the collection orders by cell, gene, umi: re-sort expected the same way
        final List<UMIKeyCount> sortedOut = new ArrayList<>(out);
        sortedOut.sort(UMIKeyCount.sortOrder(true));
        for (int i = 0; i < out.size(); i++) {
            Assert.assertSame(out.get(i), sortedOut.get(i), "output must be in comparator order");
            final UMIKeyCount k = out.get(i);
            Assert.assertEquals(k.getCount(), (int) exp.get(k.cell + "|" + k.gene + "|" + k.umi));
        }
        Assert.assertEquals(expKeys.size(), out.size());
    }

    @Test
    public void testNoSpill() {
        final List<UMIKeyCount> in = randomInput(1000, 200, 1);
        check(in, 100_000, 0);
    }

    @Test
    public void testManySpills() {
        final List<UMIKeyCount> in = randomInput(10_000, 800, 2);
        check(in, 500, 20);
    }

    @Test
    public void testAllEqualKeys() {
        final List<UMIKeyCount> in = new ArrayList<>();
        for (int i = 0; i < 5000; i++) {
            in.add(new UMIKeyCount("c", "g", "u", 1, 7 + (4999 - i)));
        }
        final List<UMIKeyCount> out = run(in, 100, 50);
        Assert.assertEquals(out.size(), 1);
        Assert.assertEquals(out.get(0).getCount(), 5000);
        Assert.assertEquals(out.get(0).getFirstSeen(), 7, "smallest firstSeen survives combining across spill files");
    }

    @Test
    public void testSpillCreatesMissingTempDir() throws IOException {
        // java.io.tmpdir is global state, so it has to be put back
        final String original = System.getProperty("java.io.tmpdir");
        final Path base = Files.createTempDirectory("combining_sort_test");
        try {
            System.setProperty("java.io.tmpdir", base.resolve("does_not_exist_yet").toString());
            final List<UMIKeyCount> in = randomInput(2000, 300, 4);
            check(in, 100, 20);
        } finally {
            System.setProperty("java.io.tmpdir", original);
        }
    }

    @Test
    public void testEmpty() {
        Assert.assertTrue(run(new ArrayList<>(), 10, 0).isEmpty());
    }

    @Test
    public void testExactMultipleOfChunkSize() {
        final List<UMIKeyCount> in = randomInput(1000, 300, 3);
        check(in, 100, 10);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testAddAfterIteratorFails() {
        final CombiningSortingCollection<UMIKeyCount> c = new CombiningSortingCollection<>(new UMIKeyCountCodec(),
                UMIKeyCount.sortOrder(true), UMIKeyCount::combine, 10);
        c.iterator();
        c.add(new UMIKeyCount("c", "g", "u", 1, 0));
    }
}
