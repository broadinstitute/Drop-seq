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
package org.broadinstitute.dropseqrna.utils.batchwriter;

import org.broadinstitute.dropseqrna.utils.ObjectSink;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Consumer;

/**
 * Collects items on the producing thread and hands them to a writer thread of its {@link BatchSinkGroup} in batches,
 * when a batch is full and when the sink is closed.  Create with {@link BatchSinkGroup#open}.
 *
 * <p>The items are kept by reference until their batch is processed, so the producer must not change or reuse an
 * item after adding it.
 *
 * @param <T> The type of the items.
 */
public final class BatchSink<T> implements ObjectSink<T> {

    private final BatchSinkGroup group;
    private final Consumer<List<T>> batchWriter;
    private final Runnable onClose;
    private final int batchSize;
    private final int maxBatchesInFlight;

    // Used by the producing thread only.
    private List<T> current;
    private boolean closed = false;

    // Guarded by lock.  inFlight counts the batches that are queued or being processed.  scheduled is true while a
    // processing task for this sink is waiting for a writer thread or running, and then the queue is not empty.
    private final Object lock = new Object();
    private final ArrayDeque<List<T>> queue = new ArrayDeque<>();
    private int inFlight = 0;
    private boolean scheduled = false;

    BatchSink(final BatchSinkGroup group, final Consumer<List<T>> batchWriter, final Runnable onClose,
              final int batchSize, final int maxBatchesInFlight) {
        this.group = group;
        this.batchWriter = batchWriter;
        this.onClose = onClose;
        this.batchSize = batchSize;
        this.maxBatchesInFlight = maxBatchesInFlight;
        this.current = new ArrayList<>(batchSize);
    }

    /**
     * Adds an item.  When this completes a batch the batch is handed over, which waits if the sink already has
     * maxBatchesInFlight batches in flight.
     * @throws IllegalStateException if the sink is closed.
     * @throws RuntimeException if processing an earlier batch failed.
     */
    @Override
    public void add(final T item) {
        if (closed) {
            throw new IllegalStateException("The sink is closed");
        }
        current.add(item);
        if (current.size() >= batchSize) {
            handOver();
        }
    }

    /**
     * Hands over the items collected so far, waits until every batch of this sink has been processed, then runs the
     * onClose action.  Does nothing if already closed.
     * @throws RuntimeException if processing any batch failed.
     */
    public void close() {
        if (closed) {
            return;
        }
        closed = true;
        if (!current.isEmpty()) {
            handOver();
        }
        synchronized (lock) {
            while (inFlight > 0) {
                await();
            }
        }
        group.throwIfFailed();
        onClose.run();
    }

    private void handOver() {
        final List<T> batch = current;
        current = new ArrayList<>(batchSize);
        synchronized (lock) {
            // After a failure the batches are discarded, so this does not wait for long.
            while (inFlight >= maxBatchesInFlight && !group.hasFailed()) {
                await();
            }
            group.throwIfFailed();
            queue.add(batch);
            inFlight++;
            if (!scheduled) {
                scheduled = true;
                group.execute(this::processOneBatch);
            }
        }
    }

    /**
     * Runs on a writer thread.  Processes one batch and then, if there are more, queues itself again at the back of
     * the pool's queue, so that one busy sink cannot keep a writer thread from the others.
     */
    private void processOneBatch() {
        final List<T> batch;
        synchronized (lock) {
            batch = queue.poll();
        }
        if (!group.hasFailed()) {
            // The failure crosses threads: it is rethrown on the producing thread, which would otherwise wait for ever.
            try {
                batchWriter.accept(batch);
            } catch (Throwable t) {
                group.fail(t);
            }
        }
        synchronized (lock) {
            inFlight--;
            if (queue.isEmpty()) {
                scheduled = false;
            } else {
                group.execute(this::processOneBatch);
            }
            lock.notifyAll();
        }
    }

    // Must be called holding lock.
    private void await() {
        try {
            lock.wait();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Interrupted while waiting for a batch to be processed", e);
        }
    }
}
