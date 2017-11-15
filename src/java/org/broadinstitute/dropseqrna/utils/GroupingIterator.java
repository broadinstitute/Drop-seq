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

import htsjdk.samtools.util.IterableOnceIterator;
import htsjdk.samtools.util.PeekableIterator;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

/**
 * Group a stream of records according to equality as defined by the given comparator.
 * It is an error in the input stream of records is not in the order defined by the comparator.
 * Each element return is a group containing all the input elements that are equal according to the comparator.
 */
public class GroupingIterator<T>  extends IterableOnceIterator<List<T>> {

    private final PeekableIterator<T> underlyingIterator;
    private final Comparator<T> comparator;

    public GroupingIterator(Iterator<T> underlyingIterator, Comparator<T> comparator) {
        this.underlyingIterator = new PeekableIterator<>(underlyingIterator);
        this.comparator = comparator;
    }

    @Override
    public boolean hasNext() {
        return underlyingIterator.hasNext();
    }

    @Override
    public List<T> next() {
        final ArrayList<T> ret = new ArrayList<>(1);
        T last = underlyingIterator.next();
        ret.add(last);
        while (underlyingIterator.hasNext()) {
            final int cmp = comparator.compare(last, underlyingIterator.peek());
            if (cmp == 0) {
                last = underlyingIterator.next();
                ret.add(last);
            } else if (cmp < 0) {
                return ret;
            } else {
                throw new IllegalStateException(String.format("Out of order iterator: %s > %s", last.toString(), underlyingIterator.peek().toString()));
            }
        }
        return ret;
    }

    @Override
    public void close() throws IOException {
        underlyingIterator.close();
        super.close();
    }
}
