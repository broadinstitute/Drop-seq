/*
 * MIT License
 *
 * Copyright 2019 Broad Institute
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

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IterableOnceIterator;

import java.util.Comparator;
import java.util.Iterator;

/**
 * Wrapper over an iterator over some type T that takes a comparator between two objects of type T
 * and will assert that the previous object is less than or equal to the next object.
 * @param <T> type of objects iterated over by underlying iterator
 * @author dmeyer
 */
public abstract class OrderAssertingIterator<T> extends IterableOnceIterator<T> {
    protected final Iterator<T> underlyingIterator;
    private final Comparator<T> orderAssertingComparator;
    private T prev = null;

    public OrderAssertingIterator(final Iterator <T> underlyingIterator,
                                  Comparator<T> orderAssertingComparator) {
        this.underlyingIterator=underlyingIterator;
        this.orderAssertingComparator=orderAssertingComparator;
    }

    @Override
    public boolean hasNext() {
        return this.underlyingIterator.hasNext();
    }

    @Override
    public T next() {
        if (this.prev == null) {
            this.prev = this.next();
            return this.prev;
        }
        T nextObj = this.underlyingIterator.next();

        if (orderAssertingComparator.compare(this.prev, nextObj) > 0)
            throw new IllegalStateException("Underlying iterator is not in asserted order!");
        return nextObj;
    }

    @Override
    public void remove() {
        this.underlyingIterator.remove();
    }

    @Override
    public void close() {
        CloserUtil.close(this.underlyingIterator);
    }
}
