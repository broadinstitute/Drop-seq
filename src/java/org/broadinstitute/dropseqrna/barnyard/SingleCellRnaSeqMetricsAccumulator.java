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
package org.broadinstitute.dropseqrna.barnyard;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.atomic.AtomicReference;

import org.broadinstitute.dropseqrna.TranscriptomeException;

import htsjdk.samtools.SAMRecord;
import picard.metrics.PerUnitMetricCollector;

/**
 * Hands each read to the metrics collector of its own cell, on the calling thread or on worker threads.
 * <p>
 * With worker threads, every cell is assigned to one worker for the whole run, so a cell's collector is only ever
 * touched by that worker: no locking and no merging of partial results.  Records are handed to workers in batches.
 * The result does not depend on the order reads arrive, so this is only for collectors whose per-cell state is
 * order independent.
 * <p>
 * Usage: call {@link #accept} for each read, then {@link #finish()} once.  {@link #finish()} must be called before
 * the collectors are read.
 */
class SingleCellRnaSeqMetricsAccumulator {

    // Batch size and queue depth made no difference in a sweep (queue 4/16/64, batch 250/1000) on the slice.
    private static final int BATCH_SIZE = 1000;
    private static final int BATCHES_IN_FLIGHT_PER_WORKER = 4;

    /** Reads handed to a worker together with the collector each one goes to (parallel arrays, no per-read objects). */
    private static class Batch {
        final PerUnitMetricCollector<?, ?, SAMRecord>[] collectors;
        final SAMRecord[] records;
        int size = 0;

        @SuppressWarnings("unchecked")
        Batch(final int capacity) {
            collectors = new PerUnitMetricCollector[capacity];
            records = new SAMRecord[capacity];
        }

        void add(final PerUnitMetricCollector<?, ?, SAMRecord> collector, final SAMRecord rec) {
            collectors[size] = collector;
            records[size++] = rec;
        }
    }

    private static final Batch END = new Batch(0);

    /** One entry per cell: the collector and the worker that owns it. */
    private static class Target {
        final PerUnitMetricCollector<?, ?, SAMRecord> collector;
        final int worker;

        Target(final PerUnitMetricCollector<?, ?, SAMRecord> collector, final int worker) {
            this.collector = collector;
            this.worker = worker;
        }
    }

    private final Map<String, Target> targets;
    private final int numWorkers;
    private final List<BlockingQueue<Batch>> queues = new ArrayList<>();
    private final List<Batch> pending = new ArrayList<>();
    private final List<Thread> threads = new ArrayList<>();
    private final AtomicReference<Throwable> failure = new AtomicReference<>();

    /**
     * @param collectorsByCell the collector of each cell to accept.  Reads of other cells are rejected by {@link #accept}.
     * @param numWorkers number of worker threads.  0 means each read is processed inline on the calling thread.
     */
    SingleCellRnaSeqMetricsAccumulator(final Map<String, ? extends PerUnitMetricCollector<?, ?, SAMRecord>> collectorsByCell,
                                       final int numWorkers) {
        this.numWorkers = numWorkers;
        this.targets = new HashMap<>(collectorsByCell.size() * 2);
        for (final Map.Entry<String, ? extends PerUnitMetricCollector<?, ?, SAMRecord>> e : collectorsByCell.entrySet()) {
            final int worker = numWorkers == 0 ? 0 : (e.getKey().hashCode() & Integer.MAX_VALUE) % numWorkers;
            targets.put(e.getKey(), new Target(e.getValue(), worker));
        }
        for (int i = 0; i < numWorkers; i++) {
            final BlockingQueue<Batch> queue = new ArrayBlockingQueue<>(BATCHES_IN_FLIGHT_PER_WORKER);
            queues.add(queue);
            pending.add(new Batch(BATCH_SIZE));
            final Thread t = new Thread(() -> work(queue), "SingleCellRnaSeqMetricsAccumulator-worker-" + i);
            t.setDaemon(true);
            threads.add(t);
            t.start();
        }
    }

    /**
     * @return true if the cell is one of the requested cells and the read was taken, false if it was ignored.
     *         The read must not be modified after this call.
     */
    boolean accept(final String cellBarcode, final SAMRecord rec) {
        final Target target = targets.get(cellBarcode);
        if (target == null)
            return false;
        if (numWorkers == 0) {
            target.collector.acceptRecord(rec);
            return true;
        }
        final Batch batch = pending.get(target.worker);
        batch.add(target.collector, rec);
        if (batch.size == BATCH_SIZE) {
            checkFailure();
            put(target.worker, batch);
            pending.set(target.worker, new Batch(BATCH_SIZE));
        }
        return true;
    }

    /** Flush, wait for all workers, and rethrow any worker failure. */
    void finish() {
        for (int i = 0; i < numWorkers; i++) {
            put(i, pending.get(i));
            put(i, END);
        }
        for (final Thread t : threads) {
            try {
                t.join();
            } catch (InterruptedException e) {
                Thread.currentThread().interrupt();
                throw new TranscriptomeException("Interrupted", e);
            }
        }
        checkFailure();
    }

    private void work(final BlockingQueue<Batch> queue) {
        try {
            while (true) {
                final Batch batch = queue.take();
                if (batch == END)
                    return;
                // After a failure keep draining so the producer never blocks.
                if (failure.get() != null)
                    continue;
                try {
                    for (int i = 0; i < batch.size; i++)
                        batch.collectors[i].acceptRecord(batch.records[i]);
                } catch (Throwable t) {
                    failure.compareAndSet(null, t);
                }
            }
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
    }

    private void checkFailure() {
        final Throwable t = failure.get();
        if (t != null)
            throw new TranscriptomeException("Failure accumulating metrics", t);
    }

    private void put(final int worker, final Batch batch) {
        try {
            queues.get(worker).put(batch);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new TranscriptomeException("Interrupted", e);
        }
    }
}
