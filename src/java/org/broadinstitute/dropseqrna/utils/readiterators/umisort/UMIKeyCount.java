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

import java.util.Comparator;

/**
 * One (cell, gene, UMI) observation with the number of reads that produced it.
 * This is what {@link org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator} sorts when it does not need to retain reads.
 * <p>
 * The count is mutable so that {@link #combine} can merge duplicates without allocating.
 * <p>
 * firstSeen is the position in the input stream of the first read that this key came from.  Sorting loses
 * the input order, and UMICollection's HashMap iteration order depends on insertion order, so
 * {@link UMIIterator} inserts the UMIs of a (cell, gene) in firstSeen order.  That gives exactly the same
 * UMICollections as sorting whole reads, which is stable, and so keeps all downstream output identical.
 */
public final class UMIKeyCount {
    public final String cell;
    public final String gene;
    public final String umi;
    private int count;
    private long firstSeen;

    public UMIKeyCount(final String cell, final String gene, final String umi, final int count, final long firstSeen) {
        this.cell = cell;
        this.gene = gene;
        this.umi = umi;
        this.count = count;
        this.firstSeen = firstSeen;
    }

    public int getCount() {
        return count;
    }

    public long getFirstSeen() {
        return firstSeen;
    }

    /** Adds the count of b (which must have equal cell, gene, UMI) to a, keeps the earlier firstSeen, and returns a. */
    public static UMIKeyCount combine(final UMIKeyCount a, final UMIKeyCount b) {
        a.count = Math.addExact(a.count, b.count);
        a.firstSeen = Math.min(a.firstSeen, b.firstSeen);
        return a;
    }

    /** Orders by the position of the first read, i.e. the order in which UMIs first appeared in the input. */
    public static Comparator<UMIKeyCount> firstSeenOrder() {
        return Comparator.comparingLong(UMIKeyCount::getFirstSeen);
    }

    /** Orders by cell then gene (cellFirst) or gene then cell, then UMI. */
    public static Comparator<UMIKeyCount> sortOrder(final boolean cellFirst) {
        return cellFirst
                ? Comparator.comparing((UMIKeyCount k) -> k.cell).thenComparing(k -> k.gene).thenComparing(k -> k.umi)
                : Comparator.comparing((UMIKeyCount k) -> k.gene).thenComparing(k -> k.cell).thenComparing(k -> k.umi);
    }

    /** Orders by cell then gene (cellFirst) or gene then cell: equal means the same UMICollection. */
    public static Comparator<UMIKeyCount> groupOrder(final boolean cellFirst) {
        return cellFirst
                ? Comparator.comparing((UMIKeyCount k) -> k.cell).thenComparing(k -> k.gene)
                : Comparator.comparing((UMIKeyCount k) -> k.gene).thenComparing(k -> k.cell);
    }

    @Override
    public String toString() {
        return "UMIKeyCount{cell=" + cell + ", gene=" + gene + ", umi=" + umi + ", count=" + count + ", firstSeen=" + firstSeen + "}";
    }
}
