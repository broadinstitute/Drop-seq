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

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.Consumer;

/**
 * A set of {@link BatchSink}s, for example one per output file, whose batches are processed by a shared pool of
 * writer threads.
 *
 * <p>The producing thread calls {@link BatchSink#add} on the sinks.  A sink collects items on the producing
 * thread, with no locking, and hands them over a whole batch at a time, so the cost of passing work between threads
 * is paid once per batch instead of once per item.  The pool has a fixed number of threads however many sinks there
 * are, so the number of threads is independent of the number of outputs.
 *
 * <p>Guarantees:
 * <ul>
 *     <li>The batches of one sink are processed in the order they were filled, and never at the same time, so a
 *     consumer that is not thread safe, such as a {@code SAMFileWriter}, can be used without synchronization.
 *     Different sinks are processed concurrently.</li>
 *     <li>Memory is bounded: a sink holds at most {@code maxBatchesInFlight} batches that have been handed over and
 *     not yet processed, plus the batch being filled.  The producer waits when that limit is reached.</li>
 *     <li>If processing a batch fails, the first failure is rethrown, wrapped in a {@link RuntimeException}, by the
 *     producer's next call that hands over a batch, or by {@link BatchSink#close} or {@link #close}.  Nothing is
 *     processed after a failure and the producer does not hang.</li>
 * </ul>
 *
 * <p>The group and its sinks must be used from one producing thread.  The writer threads are daemon threads and
 * stop by themselves if the group is never closed.
 */
public class BatchSinkGroup {

    public static final int DEFAULT_BATCH_SIZE = 1000;
    public static final int DEFAULT_MAX_BATCHES_IN_FLIGHT = 4;

    private static final AtomicInteger GROUP_COUNTER = new AtomicInteger();

    private final ThreadPoolExecutor pool;
    private final int batchSize;
    private final int maxBatchesInFlight;
    private final AtomicReference<Throwable> failure = new AtomicReference<>();
    private final List<BatchSink<?>> sinks = new ArrayList<>();

    public BatchSinkGroup(final int numWriterThreads) {
        this(numWriterThreads, DEFAULT_BATCH_SIZE, DEFAULT_MAX_BATCHES_IN_FLIGHT);
    }

    /**
     * @param numWriterThreads   Number of threads that process batches, shared by all sinks.  At least 1.
     * @param batchSize          Number of items a sink collects before it hands them over.  At least 1.
     * @param maxBatchesInFlight Per sink, the most batches that can be waiting to be processed or being processed.
     *                           At least 1.
     */
    public BatchSinkGroup(final int numWriterThreads, final int batchSize, final int maxBatchesInFlight) {
        if (numWriterThreads < 1 || batchSize < 1 || maxBatchesInFlight < 1) {
            throw new IllegalArgumentException(String.format(
                    "numWriterThreads, batchSize and maxBatchesInFlight must all be at least 1: %d, %d, %d",
                    numWriterThreads, batchSize, maxBatchesInFlight));
        }
        this.batchSize = batchSize;
        this.maxBatchesInFlight = maxBatchesInFlight;
        final int groupId = GROUP_COUNTER.incrementAndGet();
        final AtomicInteger threadCounter = new AtomicInteger();
        this.pool = new ThreadPoolExecutor(numWriterThreads, numWriterThreads, 1, TimeUnit.SECONDS,
                new LinkedBlockingQueue<>(), runnable -> {
            final Thread thread = new Thread(runnable, "BatchSinkWriter-" + groupId + "-" + threadCounter.incrementAndGet());
            thread.setDaemon(true);
            return thread;
        });
        // Idle threads of a group that was not closed because of a failure go away by themselves.
        this.pool.allowCoreThreadTimeOut(true);
    }

    /**
     * @param batchWriter Processes one batch.  Called on a writer thread, for the batches of this sink in order and
     *                    never concurrently with itself.
     * @param onClose     Run on the thread that calls {@link BatchSink#close} once every batch has been processed,
     *                    for example to close the file being written.  Not run if a batch failed.
     */
    public <T> BatchSink<T> open(final Consumer<List<T>> batchWriter, final Runnable onClose) {
        final BatchSink<T> sink = new BatchSink<>(this, batchWriter, onClose, batchSize, maxBatchesInFlight);
        sinks.add(sink);
        return sink;
    }

    /**
     * Closes every sink that is still open, in the order they were opened, then stops the writer threads.
     * @throws RuntimeException if processing any batch failed.
     */
    public void close() {
        for (final BatchSink<?> sink : sinks) {
            sink.close();
        }
        pool.shutdown();
    }

    void execute(final Runnable task) {
        pool.execute(task);
    }

    boolean hasFailed() {
        return failure.get() != null;
    }

    void fail(final Throwable t) {
        failure.compareAndSet(null, t);
    }

    void throwIfFailed() {
        final Throwable t = failure.get();
        if (t != null) {
            throw new RuntimeException("Failed processing a batch: " + t, t);
        }
    }
}
