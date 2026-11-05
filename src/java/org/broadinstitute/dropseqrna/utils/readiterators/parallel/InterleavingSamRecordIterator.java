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
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * Merges several unsorted inputs by taking up to a fixed number of records from each in turn, round robin.
 * <p>
 * Unlike {@link org.broadinstitute.dropseqrna.utils.readiterators.UnsortedMergingSamRecordIterator}, which does not
 * ask the next reader for an iterator until the previous reader is exhausted, this asks every reader for its iterator
 * up front.  For readers created by a {@link BamDecodePool} that starts reading every file at once, so all of the
 * pool's threads can be used however few files remain, instead of one file at a time.
 * <p>
 * The order is deterministic (it depends only on the chunk size and the files), but it is not the order of the files,
 * so use this only where record order does not matter.  A failure reading one input is thrown when that input's turn
 * comes, after the records of that input that came before the failure.
 */
class InterleavingSamRecordIterator implements CloseableIterator<SAMRecord> {
    private final SAMFileHeader outputHeader;
    private final List<SamReader> readers;
    private final int chunkSize;
    /** Iterators of the inputs that may still have records, in rotation order. */
    private final List<Iterator<SAMRecord>> active;
    /** Index in active of the input whose turn it is. */
    private int current = 0;
    /** Records taken from the current input in this turn. */
    private int emitted = 0;

    /**
     * @param outputHeader If non-null, each SAMRecord has its header set to this.
     * @param chunkSize records taken from one input before moving to the next, at least 1.  Using the pool's batch
     *                  size lines the turns up with the batches the pool decodes.
     */
    InterleavingSamRecordIterator(final SAMFileHeader outputHeader, final List<SamReader> readers, final int chunkSize) {
        if (chunkSize < 1) {
            throw new IllegalArgumentException("chunkSize must be at least 1: " + chunkSize);
        }
        this.outputHeader = outputHeader;
        this.readers = readers;
        this.chunkSize = chunkSize;
        this.active = new ArrayList<>(readers.size());
        try {
            for (final SamReader reader : readers) {
                active.add(reader.iterator());
            }
        } catch (RuntimeException | Error e) {
            // Don't leave the files that already started reading holding the pool busy.
            CloserUtil.close(readers);
            throw e;
        }
    }

    @Override
    public boolean hasNext() {
        while (!active.isEmpty()) {
            if (emitted >= chunkSize) {
                current = (current + 1) % active.size();
                emitted = 0;
            }
            if (active.get(current).hasNext()) {
                return true;
            }
            // This input is finished.  The next one slides into its index, so current stays where it is.
            active.remove(current);
            if (current >= active.size()) {
                current = 0;
            }
            emitted = 0;
        }
        return false;
    }

    @Override
    public SAMRecord next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        final SAMRecord rec = active.get(current).next();
        emitted++;
        if (outputHeader != null) {
            rec.setHeader(outputHeader);
        }
        return rec;
    }

    @Override
    public void close() {
        CloserUtil.close(readers);
    }
}
