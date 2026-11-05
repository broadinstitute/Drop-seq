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
package org.broadinstitute.dropseqrna.utils.readiterators.parallel;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import java.io.Closeable;
import java.nio.file.Path;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Semaphore;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.IntConsumer;

/**
 * A fixed pool of threads that reads and decodes BAM records for any number of files.
 * <p>
 * The readers returned by {@link #open(Path)} have no threads of their own: all reading and decoding is done by
 * tasks on this pool, so the thread count stays at {@code numThreads} however many files are open.  Nothing is
 * read, and no thread is started, until {@link SamReader#iterator()} is called on a reader.
 * <p>
 * For when more threads help, and when they do not, see the numThreads documentation of
 * {@link org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil}.
 * <p>
 * Only plain BAM files are read in parallel.  Anything else (SAM, CRAM, ...) is opened with the
 * {@link SamReaderFactory} given to the pool.  Records are always fully decoded on the pool threads, as with
 * {@link SamReaderFactory.Option#EAGERLY_DECODE}, and validated with the factory's validation stringency.  The
 * factory's other options are not applied to the BAMs read in parallel.
 */
public class BamDecodePool implements Closeable {
    /**
     * Records per batch.  Together with the read-ahead (by default numThreads batches, shared by all files) this keeps
     * a small, fixed number of records in flight: measured throughput was best with about 1,000 to 4,000 records in
     * flight, fell off with 100-record batches (task overhead) and with 8,000 or more (a cliff, probably cache and
     * garbage collection pressure).
     */
    static final int DEFAULT_BATCH_SIZE = 250;
    private static final AtomicInteger POOL_COUNTER = new AtomicInteger();

    private final SamReaderFactory factory;
    private final SamReaderFactory sourceFactory =
            SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
    private final ExecutorService executor;
    private final int batchSize;
    private final Semaphore readAheadPermits;
    private final IntConsumer beforeBatchDecode;
    private final Set<ParallelBamReader> openReaders = ConcurrentHashMap.newKeySet();
    private volatile boolean closed = false;

    public BamDecodePool(final int numThreads) {
        this(numThreads, SamReaderFactory.makeDefault());
    }

    /**
     * @param numThreads number of pool threads, at least 1
     * @param factory used to read headers, and to open anything that is not a plain BAM
     */
    public BamDecodePool(final int numThreads, final SamReaderFactory factory) {
        this(numThreads, factory, DEFAULT_BATCH_SIZE, numThreads, null);
    }

    /**
     * @param batchSize records per batch
     * @param readAhead batches that may be read and decoded ahead of the consumer, beyond the one batch every open
     *                  file always has outstanding.  Shared by all the files of the pool, so memory is bounded by about
     *                  (files + readAhead) batches however many files are open.
     * @param beforeBatchDecode test hook, called on the pool thread with the batch's index in its file just before
     *                          the batch is decoded.  May be null.
     */
    BamDecodePool(final int numThreads, final SamReaderFactory factory, final int batchSize, final int readAhead,
                  final IntConsumer beforeBatchDecode) {
        if (numThreads < 1) {
            throw new IllegalArgumentException("numThreads must be at least 1: " + numThreads);
        }
        if (batchSize < 1 || readAhead < 0) {
            throw new IllegalArgumentException("batchSize must be at least 1 and readAhead at least 0");
        }
        this.factory = factory;
        this.batchSize = batchSize;
        this.readAheadPermits = new Semaphore(readAhead);
        this.beforeBatchDecode = beforeBatchDecode;
        final int poolId = POOL_COUNTER.incrementAndGet();
        final AtomicInteger workerCounter = new AtomicInteger();
        final ThreadFactory threadFactory = r -> {
            final Thread t = new Thread(r, "BamDecodePool-" + poolId + "-worker-" + workerCounter.incrementAndGet());
            t.setDaemon(true);
            return t;
        };
        // Core size == max size: threads are created as tasks arrive, so an unused pool starts none.
        this.executor = Executors.newFixedThreadPool(numThreads, threadFactory);
    }

    /**
     * Open a file.  A plain BAM gets a reader that decodes on this pool; anything else is returned as opened by the
     * pool's SamReaderFactory.  Only the header is read here, and no file stays open until the reader's iterator is
     * requested.
     */
    public SamReader open(final Path path) {
        if (closed) {
            throw new IllegalStateException("Pool is closed");
        }
        final SamReader opened = factory.open(path);
        if (opened.type() != SamReader.Type.BAM_TYPE) {
            return opened;
        }
        final SAMFileHeader header = opened.getFileHeader();
        CloserUtil.close(opened);
        final ParallelBamReader reader = new ParallelBamReader(this, path, header, factory.validationStringency());
        openReaders.add(reader);
        return reader;
    }

    /** Closes every reader this pool opened and stops the threads. */
    @Override
    public void close() {
        closed = true;
        for (final ParallelBamReader reader : openReaders) {
            reader.close();
        }
        executor.shutdownNow();
        try {
            executor.awaitTermination(30, TimeUnit.SECONDS);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
        }
    }

    /**
     * Merge unsorted readers opened by this pool so that all of them are read at the same time, taking one batch from
     * each in turn.  Records are not in the order of the files; use only where order does not matter.  Closing the
     * returned iterator closes the readers.
     * @param outputHeader If non-null, each record has its header set to this.
     */
    public CloseableIterator<SAMRecord> interleave(final SAMFileHeader outputHeader, final List<SamReader> readers) {
        return new InterleavingSamRecordIterator(outputHeader, readers, batchSize);
    }

    /** @return an iterator over the same records that also closes this pool when it is closed. */
    public CloseableIterator<SAMRecord> closeWhenClosed(final CloseableIterator<SAMRecord> iterator) {
        return new CloseableIterator<>() {
            @Override
            public boolean hasNext() {
                return iterator.hasNext();
            }

            @Override
            public SAMRecord next() {
                return iterator.next();
            }

            @Override
            public void close() {
                iterator.close();
                BamDecodePool.this.close();
            }
        };
    }

    ExecutorService executor() {
        return executor;
    }

    int batchSize() {
        return batchSize;
    }

    /** @return true if the caller may have one more batch outstanding beyond its guaranteed one; release it later. */
    boolean tryAcquireReadAhead() {
        return readAheadPermits.tryAcquire();
    }

    void releaseReadAhead(final int permits) {
        if (permits > 0) {
            readAheadPermits.release(permits);
        }
    }

    /** Opens the private htsjdk reader of each BAM: lazy record decoding, and validation left to the pool. */
    SamReaderFactory sourceFactory() {
        return sourceFactory;
    }

    IntConsumer beforeBatchDecode() {
        return beforeBatchDecode;
    }

    void readerClosed(final ParallelBamReader reader) {
        openReaders.remove(reader);
    }
}
