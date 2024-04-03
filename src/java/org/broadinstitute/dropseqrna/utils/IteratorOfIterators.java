/*
 * MIT License
 *
 * Copyright 2023 Broad Institute
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

import htsjdk.samtools.util.CloseableIterator;

import java.util.*;
import java.util.function.Function;
import java.util.stream.StreamSupport;

/**
 * Converts an iterator of iterators of T into an iterator of T
 */
public class IteratorOfIterators<T>
        implements CloseableIterator<T> {
    private final Iterator<Iterator<T>> iteratorIt;
    private Iterator<T> currentIt;

    public IteratorOfIterators(Iterator<Iterator<T>> iteratorIt) {
        this.iteratorIt = iteratorIt;
        if (this.iteratorIt.hasNext()) {
            currentIt = this.iteratorIt.next();
        } else {
            currentIt = Collections.emptyIterator();
        }
    }

    @Override
    public boolean hasNext() {
        while (!currentIt.hasNext()) {
            if (currentIt instanceof CloseableIterator<T>) {
                ((CloseableIterator<T>) currentIt).close();
            }

            if (iteratorIt.hasNext()) {
                currentIt = iteratorIt.next();
            } else {
                currentIt = Collections.emptyIterator();
                return false;
            }
        }
        return true;
    }

    @Override
    public T next() {
        if (!hasNext()) {
            throw new IllegalStateException("Cannot call next() when !hasNext()");
        }
        return currentIt.next();
    }

    @Override
    public void close() {
        if (currentIt instanceof CloseableIterator<T>) {
            ((CloseableIterator<T>) currentIt).close();
        }
        while (iteratorIt.hasNext()) {
            final Iterator<T> it = iteratorIt.next();
            if (it instanceof CloseableIterator<T>) {
                ((CloseableIterator<T>) it).close();
            }
        }
    }

    public static <T, S> IteratorOfIterators<T> fromIterator(Iterator<S> outerIter, Function<S, Iterator<T>> mapIter) {
        final Iterator<Iterator<T>> iteratorIt =
                StreamSupport
                        .stream(Spliterators.spliteratorUnknownSize(outerIter, Spliterator.ORDERED), false)
                        .map(mapIter)
                        .iterator();
        return new IteratorOfIterators<>(iteratorIt);
    }
}
