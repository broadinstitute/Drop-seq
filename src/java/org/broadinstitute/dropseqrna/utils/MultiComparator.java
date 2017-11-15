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
package org.broadinstitute.dropseqrna.utils;

import java.util.Comparator;
import java.util.List;

/**
 * Create a comparator from a list of comparators that all operate on the same type.  The first comparator is the
 * primary sort, second comparator is secondary sort, etc.
 */
public class MultiComparator<T> implements Comparator<T> {

    private final Comparator<T>[] comparators;

    @SuppressWarnings("unchecked")
    public MultiComparator(final List<Comparator<T>> comparators) {
        this.comparators = comparators.toArray(new Comparator[comparators.size()]);
    }

    public MultiComparator(final Comparator<T>...comparators) {
        this.comparators = comparators;
    }

    @Override
    public int compare(final T o1, final T o2) {
        for (final Comparator<T> comparator : comparators) {
            final int cmp = comparator.compare(o1, o2);
            if (cmp != 0)
				return cmp;
        }
        return 0;
    }
}
