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

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IterableOnceIterator;

import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Queue;

/**
 * Iterator wrapper that allows concrete implementations to produce 0, 1 or more output records for each input
 * record. processRecord() must be implemented, and it must call queueRecordForOutput() for each output record
 * to be emitted.
 *
 * Note that the decision about what to emit can only be based on the current record passed to processRecord.
 */
public abstract class CountChangingIteratorWrapper<T> extends IterableOnceIterator<T> {

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
    public void close() throws IOException {
        CloserUtil.close(underlyingIterator);
        super.close();
    }
}
