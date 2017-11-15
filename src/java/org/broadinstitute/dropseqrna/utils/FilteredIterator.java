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

import java.util.Iterator;

import htsjdk.samtools.util.IterableOnceIterator;
import htsjdk.samtools.util.PeekableIterator;

public abstract class FilteredIterator<T> extends IterableOnceIterator<T> {

    final PeekableIterator<T> it;
    boolean firstTime = true;

    protected FilteredIterator(final Iterator<T> underlyingIterator) {
        it = new PeekableIterator<>(underlyingIterator);
    }

    /**
     * @return true if record should be skipped
     */
    public abstract boolean filterOut(final T rec);

    private void skipUndesired() {
        while (it.hasNext() && filterOut(it.peek()))
			it.next();
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
