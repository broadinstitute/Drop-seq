/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2015 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
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
