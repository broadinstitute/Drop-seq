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

import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import java.nio.file.Path;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.concurrent.CancellationException;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;

/**
 * Reads BAM files with the decoding of their records spread over several threads, to speed up scanning every record.
 * One or many files can be read, and all of them share the same threads.  Each file has its own instance of this class.
 * Used by {@link org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil} when more than one thread is
 * requested.
 *
 * <h2>When more threads help</h2>
 * Based on benchmarks, so check it against your own data.
 * <ul>
 * <li><b>Several input files: reading scales with many threads.</b>  The files are read at the same time, so the gain
 * follows the number of files, roughly one thread per file.</li>
 * <li><b>One aligned file: reading scales with up to 2 threads.</b>  An uncompressed file kept gaining up to about
 * 4 threads.</li>
 * <li><b>One unaligned file: asynchronous IO is better.</b>  Such records are cheap to decode, so extra threads gain
 * little, and htsjdk's asynchronous IO is as fast or faster.  SamFileMergeUtil uses it for a single unaligned
 * input.</li>
 * </ul>
 * Only reading and decoding are done on the extra threads.  Work the caller does on each record is not sped up.
 */
class ParallelBamReader implements SamReader {
    private final BamDecodePool pool;
    private final Path path;
    private final SAMFileHeader header;
    private final ValidationStringency stringency;
    private final int batchSize;

    private final Object lock = new Object();
    // All of the following are guarded by lock.
    private final ArrayDeque<CompletableFuture<Decoded>> queue = new ArrayDeque<>();
    private int outstanding = 0;
    /** Pool read-ahead permits held: always max(0, outstanding - 1), since one batch per file is guaranteed. */
    private int extraHeld = 0;
    private boolean readInFlight = false;
    private boolean eof = false;
    private volatile boolean closed = false;   // written under lock; volatile so iterators can check it cheaply
    private boolean iteratorCreated = false;
    private int nextBatchIndex = 0;

    /** The private htsjdk reader.  Only the (single at a time) read task uses the iterator; close() closes it. */
    private SamReader source;
    private CloseableIterator<SAMRecord> sourceIterator;
    /** Records handed to decode tasks so far.  Only touched by the read task. */
    private long recordsRead = 0;

    /**
     * @param header the header to return from {@link #getFileHeader()}; the caller may change it
     * @param stringency validation stringency to apply to each record
     */
    ParallelBamReader(final BamDecodePool pool, final Path path, final SAMFileHeader header,
                      final ValidationStringency stringency) {
        this.pool = pool;
        this.path = path;
        this.header = header;
        this.stringency = stringency;
        this.batchSize = pool.batchSize();
    }

    @Override
    public SAMFileHeader getFileHeader() {
        return header;
    }

    @Override
    public Type type() {
        return Type.BAM_TYPE;
    }

    @Override
    public String getResourceDescription() {
        return path.toUri().toString();
    }

    @Override
    public boolean hasIndex() {
        return false;
    }

    @Override
    public Indexing indexing() {
        throw new UnsupportedOperationException("Not supported by a parallel BAM reader");
    }

    @Override
    public SAMRecordIterator iterator() {
        synchronized (lock) {
            if (closed) {
                throw new IllegalStateException("Reader is closed: " + path);
            }
            if (iteratorCreated) {
                throw new IllegalStateException("iterator() may only be called once: " + path);
            }
            iteratorCreated = true;
        }
        source = pool.sourceFactory().open(path);
        sourceIterator = source.iterator();
        schedule();
        return new SamReader.AssertingIterator(new BatchIterator());
    }

    @Override
    public void close() {
        final List<CompletableFuture<Decoded>> toCancel;
        synchronized (lock) {
            if (closed) {
                return;
            }
            closed = true;
            toCancel = new ArrayList<>(queue);
            queue.clear();
            outstanding = 0;
            pool.releaseReadAhead(extraHeld);
            extraHeld = 0;
        }
        toCancel.forEach(f -> f.cancel(false));
        CloserUtil.close(sourceIterator);
        CloserUtil.close(source);
        pool.readerClosed(this);
    }

    @Override
    public SAMRecordIterator query(String sequence, int start, int end, boolean contained) {
        throw unsupported();
    }

    @Override
    public SAMRecordIterator queryOverlapping(String sequence, int start, int end) {
        throw unsupported();
    }

    @Override
    public SAMRecordIterator queryContained(String sequence, int start, int end) {
        throw unsupported();
    }

    @Override
    public SAMRecordIterator query(QueryInterval[] intervals, boolean contained) {
        throw unsupported();
    }

    @Override
    public SAMRecordIterator queryOverlapping(QueryInterval[] intervals) {
        throw unsupported();
    }

    @Override
    public SAMRecordIterator queryContained(QueryInterval[] intervals) {
        throw unsupported();
    }

    @Override
    public SAMRecordIterator queryUnmapped() {
        throw unsupported();
    }

    @Override
    public SAMRecordIterator queryAlignmentStart(String sequence, int alignmentStart) {
        throw unsupported();
    }

    @Override
    public SAMRecord queryMate(SAMRecord rec) {
        throw unsupported();
    }

    private static UnsupportedOperationException unsupported() {
        return new UnsupportedOperationException("Not supported by a parallel BAM reader");
    }

    // ---------------------------------------------------------------------------------------------------------
    // Scheduling

    /** Start a read of the next batch if this file has room for one and no read is running. */
    private void schedule() {
        final CompletableFuture<Decoded> future;
        final int batchIndex;
        synchronized (lock) {
            if (closed || eof || readInFlight) {
                return;
            }
            if (outstanding > 0) {
                // The first batch is always allowed, so a merge can never be starved of the next batch it needs.
                // Further read-ahead comes from the pool-wide allowance.
                if (!pool.tryAcquireReadAhead()) {
                    return;
                }
                extraHeld++;
            }
            readInFlight = true;
            outstanding++;
            future = new CompletableFuture<>();
            queue.add(future);   // queue order is file order
            batchIndex = nextBatchIndex++;
        }
        pool.executor().execute(() -> readTask(future, batchIndex));
    }

    /** The records of one batch, as read and not yet decoded, and how the read ended. */
    private static class RawBatch {
        final List<SAMRecord> records;
        final boolean atEof;
        /** 1-based number in the file of the first record, for validation messages. */
        final long firstRecordNumber;
        /** Set if reading failed after the records in this batch. */
        final Throwable readError;

        RawBatch(List<SAMRecord> records, boolean atEof, long firstRecordNumber, Throwable readError) {
            this.records = records;
            this.atEof = atEof;
            this.firstRecordNumber = firstRecordNumber;
            this.readError = readError;
        }
    }

    private void readTask(final CompletableFuture<Decoded> future, final int batchIndex) {
        try {
            final RawBatch raw = readBatch();
            synchronized (lock) {
                readInFlight = false;
                if (raw.atEof || raw.readError != null) {
                    eof = true;
                }
            }
            schedule();   // overlap the next read with this batch's decode
            pool.executor().execute(() -> decodeTask(raw, future, batchIndex));
        } catch (Throwable t) {
            // Whatever went wrong, the consumer is waiting on this future.
            synchronized (lock) {
                readInFlight = false;
                eof = true;
            }
            future.completeExceptionally(t);
        }
    }

    /**
     * Pull up to batchSize records from htsjdk's reader.  A failure, including an Error such as OutOfMemoryError from
     * a corrupt record, is returned with the records read before it, so it reaches the consumer in order.
     */
    private RawBatch readBatch() {
        final List<SAMRecord> records = new ArrayList<>(batchSize);
        boolean atEof = false;
        Throwable readError = null;
        try {
            while (records.size() < batchSize) {
                if (!sourceIterator.hasNext()) {
                    atEof = true;
                    break;
                }
                records.add(sourceIterator.next());
            }
        } catch (Throwable t) {
            readError = t;
        }
        final long firstRecordNumber = recordsRead + 1;
        recordsRead += records.size();
        return new RawBatch(records, atEof, firstRecordNumber, readError);
    }

    /**
     * Decode and validate a batch.  A failure on a record does not discard the records before it: they are
     * delivered first and the failure is thrown after them, as with htsjdk's own reader.
     */
    private void decodeTask(final RawBatch raw, final CompletableFuture<Decoded> future, final int batchIndex) {
        if (future.isCancelled()) {
            return;
        }
        final List<SAMRecord> good = new ArrayList<>(raw.records.size());
        Throwable error = null;
        try {
            if (pool.beforeBatchDecode() != null) {
                pool.beforeBatchDecode().accept(batchIndex);
            }
            for (final SAMRecord rec : raw.records) {
                decodeFully(rec);
                if (stringency != ValidationStringency.SILENT) {
                    SAMUtils.processValidationErrors(rec.isValid(stringency == ValidationStringency.STRICT),
                            raw.firstRecordNumber + good.size(), stringency);
                }
                good.add(rec);
            }
            error = raw.readError;
        } catch (Throwable t) {
            error = t;
        }
        future.complete(new Decoded(good, error));
    }

    /** Touch everything htsjdk decodes lazily, using only public accessors (what EAGERLY_DECODE amounts to). */
    private static void decodeFully(final SAMRecord rec) {
        rec.getReadName();
        rec.getCigar();
        rec.getReadBases();
        rec.getBaseQualities();
        rec.getAttributes();
    }

    /** The result of decoding a batch: its records, then an optional failure that comes after them. */
    private static class Decoded {
        final List<SAMRecord> records;
        final Throwable error;

        Decoded(final List<SAMRecord> records, final Throwable error) {
            this.records = records;
            this.error = error;
        }
    }

    // ---------------------------------------------------------------------------------------------------------
    // Consuming

    private class BatchIterator implements CloseableIterator<SAMRecord> {
        private Iterator<SAMRecord> current = Collections.emptyIterator();
        /** Thrown once the records of the batch it came from have been handed out. */
        private Throwable pendingError = null;

        @Override
        public boolean hasNext() {
            if (closed) {
                return false;
            }
            while (!current.hasNext()) {
                if (pendingError != null) {
                    throwUnchecked(pendingError);
                }
                final CompletableFuture<Decoded> head;
                synchronized (lock) {
                    head = queue.peek();
                }
                if (head == null) {
                    return false;   // schedule() always runs after a poll, so an empty queue means end of file (or closed)
                }
                final Decoded batch;
                try {
                    batch = head.get();
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                    throw new RuntimeException("Interrupted while reading " + path, e);
                } catch (CancellationException e) {
                    return false;   // closed
                } catch (ExecutionException e) {
                    throwUnchecked(e.getCause());
                    return false;   // not reached
                }
                synchronized (lock) {
                    if (queue.peek() == head) {
                        queue.poll();
                        outstanding--;
                        final int shouldHold = Math.max(0, outstanding - 1);
                        pool.releaseReadAhead(extraHeld - shouldHold);
                        extraHeld = shouldHold;
                    }
                }
                schedule();
                current = batch.records.iterator();
                pendingError = batch.error;
            }
            return true;
        }

        private void throwUnchecked(final Throwable t) {
            if (t instanceof RuntimeException re) {
                throw re;
            }
            if (t instanceof Error err) {
                throw err;
            }
            throw new RuntimeException("Error reading " + path, t);
        }

        @Override
        public SAMRecord next() {
            if (!hasNext()) {
                throw new NoSuchElementException();
            }
            final SAMRecord rec = current.next();
            rec.setHeader(header);
            return rec;
        }

        @Override
        public void close() {
            ParallelBamReader.this.close();
        }
    }
}
