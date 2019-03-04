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

import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Queue;

import org.broadinstitute.dropseqrna.TranscriptomeException;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IterableOnceIterator;

/**
 * Iterator wrapper that allows concrete implementations to produce 0, 1 or more output records for each input
 * record. processRecord() must be implemented, and it must call queueRecordForOutput() for each output record
 * to be emitted.
 *
 * Note that the decision about what to emit can only be based on the current record passed to processRecord.
 */
public abstract class CountChangingIteratorWrapper<T> extends IterableOnceIterator<T> implements CloseableIterator<T>{

    private final Iterator<T> underlyingIterator;
    private final Queue<T> outputQueue = new LinkedList<>();

    protected CountChangingIteratorWrapper(final Iterator<T> underlyingIterator) {
        this.underlyingIterator = underlyingIterator;
        // Cannot call populateOutputQueue here, because derived object has not be constructed yet.
    }

    private void populateOutputQueue() {
        while (outputQueue.isEmpty() && underlyingIterator.hasNext()) {
            processRecord(underlyingIterator.next());
        }
    }

    protected abstract void processRecord(final T rec);

    /**
     * Called by concrete implementation to put a record on the output queue
     */
    protected void queueRecordForOutput(final T rec) {
        outputQueue.add(rec);
    }

    @Override
    public boolean hasNext() {
        populateOutputQueue();
        return !outputQueue.isEmpty();
    }

    @Override
    public T next() {
        // Called here just in case next() is called without hasNext()
        populateOutputQueue();
        return outputQueue.remove();
    }

    @Override
    public void close() {
        CloserUtil.close(underlyingIterator);
        try {
			super.close();
		} catch (IOException e) {
			throw new TranscriptomeException(e.getMessage(), e.getCause());
		}
    }
}
