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

import java.util.Iterator;

public abstract class FilteredIterator<T> extends IterableOnceIterator<T> {

    final PeekableIterator<T> it;
    boolean firstTime = true;

    protected FilteredIterator(final Iterator<T> underlyingIterator) {
        it = new PeekableIterator<>(underlyingIterator);
    }

    /**
     * @return true if record should be skipped
     */
    protected abstract boolean filterOut(final T rec);

    private void skipUndesired() {
        while (it.hasNext() && filterOut(it.peek())) {
            it.next();
        }
    }

    private void maybeSkipFirstTime() {
        if (firstTime) {
            skipUndesired();
            firstTime = false;
        }
    }

    @Override
    public void close() {
        it.close();
    }

    @Override
    public boolean hasNext() {
        maybeSkipFirstTime();
        return it.hasNext();
    }

    @Override
    public T next() {
        maybeSkipFirstTime();
        final T ret = it.next();
        skipUndesired();
        return ret;
    }

    @Override
    public void remove() {
        it.remove();
    }
}
