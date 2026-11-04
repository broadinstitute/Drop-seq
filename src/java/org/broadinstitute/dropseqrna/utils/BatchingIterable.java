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
package org.broadinstitute.dropseqrna.utils;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Presents the items of an iterator as a sequence of lists ("batches") of up to batchSize items, in the original order.
 * Every batch has exactly batchSize items except possibly the last, which has between 1 and batchSize.  There are no
 * empty batches, so an empty source gives no batches.
 *
 * <p>Intended for handing work to threads in chunks big enough that the cost of passing each chunk
 * (a queue operation, a task submission) is small compared with the work done on it.  See OrderedConcurrentMapper.
 *
 * <p>Can only be iterated once, since it consumes the source iterator.
 *
 * @param <T> The type of the items being batched.
 */
public class BatchingIterable<T> implements Iterable<List<T>> {

    private final Iterator<T> source;
    private final int batchSize;
    private boolean iterated = false;

    /**
     * @param source    The items to batch.
     * @param batchSize The maximum number of items in a batch.  Must be at least 1.
     */
    public BatchingIterable(final Iterator<T> source, final int batchSize) {
        if (batchSize < 1) {
            throw new IllegalArgumentException("batchSize must be at least 1: " + batchSize);
        }
        this.source = source;
        this.batchSize = batchSize;
    }

    /**
     * @throws IllegalStateException if called more than once.
     */
    @Override
    public Iterator<List<T>> iterator() {
        if (this.iterated) {
            throw new IllegalStateException("A BatchingIterable can only be iterated once.");
        }
        this.iterated = true;
        return new Iterator<>() {
            @Override
            public boolean hasNext() {
                return source.hasNext();
            }

            @Override
            public List<T> next() {
                if (!source.hasNext()) {
                    throw new NoSuchElementException();
                }
                final List<T> batch = new ArrayList<>(batchSize);
                while (batch.size() < batchSize && source.hasNext()) {
                    batch.add(source.next());
                }
                return batch;
            }
        };
    }
}
