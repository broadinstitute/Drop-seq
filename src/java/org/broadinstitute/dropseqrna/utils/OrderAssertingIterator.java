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
