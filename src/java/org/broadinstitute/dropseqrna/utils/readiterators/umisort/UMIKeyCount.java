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
 */
public final class UMIKeyCount {
    private static final Comparator<UMIKeyCount> CELL_FIRST_ORDER = (a, b) -> {
        int c = a.cell.compareTo(b.cell);
        if (c != 0) return c;
        c = a.gene.compareTo(b.gene);
        return c != 0 ? c : a.umi.compareTo(b.umi);
    };
    private static final Comparator<UMIKeyCount> GENE_FIRST_ORDER = (a, b) -> {
        int c = a.gene.compareTo(b.gene);
        if (c != 0) return c;
        c = a.cell.compareTo(b.cell);
        return c != 0 ? c : a.umi.compareTo(b.umi);
    };
    private static final Comparator<UMIKeyCount> CELL_FIRST_GROUP = (a, b) -> {
        final int c = a.cell.compareTo(b.cell);
        return c != 0 ? c : a.gene.compareTo(b.gene);
    };
    private static final Comparator<UMIKeyCount> GENE_FIRST_GROUP = (a, b) -> {
        final int c = a.gene.compareTo(b.gene);
        return c != 0 ? c : a.cell.compareTo(b.cell);
    };

    public final String cell;
    public final String gene;
    public final String umi;
    private int count;

    public UMIKeyCount(final String cell, final String gene, final String umi, final int count) {
        this.cell = cell;
        this.gene = gene;
        this.umi = umi;
        this.count = count;
    }

    public int getCount() {
        return count;
    }

    /** Adds the count of b (which must have equal cell, gene, UMI) to a, and returns a. */
    public static UMIKeyCount combine(final UMIKeyCount a, final UMIKeyCount b) {
        a.count = Math.addExact(a.count, b.count);
        return a;
    }

    /** Orders by cell then gene (cellFirst) or gene then cell, then UMI. */
    public static Comparator<UMIKeyCount> sortOrder(final boolean cellFirst) {
        return cellFirst ? CELL_FIRST_ORDER : GENE_FIRST_ORDER;
    }

    /** Orders by cell then gene (cellFirst) or gene then cell: equal means the same UMICollection. */
    public static Comparator<UMIKeyCount> groupOrder(final boolean cellFirst) {
        return cellFirst ? CELL_FIRST_GROUP : GENE_FIRST_GROUP;
    }

    @Override
    public String toString() {
        return "UMIKeyCount{cell=" + cell + ", gene=" + gene + ", umi=" + umi + ", count=" + count + "}";
    }
}
