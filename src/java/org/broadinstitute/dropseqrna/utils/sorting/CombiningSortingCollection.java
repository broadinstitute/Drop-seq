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
package org.broadinstitute.dropseqrna.utils.sorting;

import htsjdk.samtools.Defaults;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.TempStreamFactory;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.PriorityQueue;
import java.util.function.BinaryOperator;

/**
 * Like htsjdk's SortingCollection, but records that compare equal are merged into one with a combiner,
 * both when a chunk is spilled and when the spilled chunks are merged. Use it when records carry a count
 * and many are duplicates, so that far less is written to disk and read back.
 * <p>
 * Usage: add() everything, then call iterator() once. The combiner may mutate and return its first argument.
 * Records that compare equal must be interchangeable apart from what the combiner merges.
 * The result is in comparator order with no two equal records.
 * <p>
 * Temp files go to the htsjdk default temp dir (set by the TMP_DIR argument of command line programs).
 * They are deleted when the iterator is exhausted or closed, and on JVM exit.
 */
public class CombiningSortingCollection<T> {

    private final SortingCollection.Codec<T> codec;
    private final Comparator<T> comparator;
    private final BinaryOperator<T> combiner;
    private final int maxRecordsInRam;
    private final TempStreamFactory tempStreamFactory = new TempStreamFactory();
    private final List<File> files = new ArrayList<>();
    private final ArrayList<T> ram;
    private boolean iterating = false;

    /**
     * @param codec writes and reads records of spilled chunks
     * @param comparator defines the order, and which records are merged (compare == 0)
     * @param combiner merges two equal records into one
     * @param maxRecordsInRam number of records (before combining) held in memory before a chunk is spilled
     */
    public CombiningSortingCollection(final SortingCollection.Codec<T> codec, final Comparator<T> comparator,
                                      final BinaryOperator<T> combiner, final int maxRecordsInRam) {
        if (maxRecordsInRam < 1) {
            throw new IllegalArgumentException("maxRecordsInRam must be positive: " + maxRecordsInRam);
        }
        this.codec = codec;
        this.comparator = comparator;
        this.combiner = combiner;
        this.maxRecordsInRam = maxRecordsInRam;
        this.ram = new ArrayList<>(Math.min(maxRecordsInRam, 1 << 16));
    }

    public void add(final T record) {
        if (iterating) {
            throw new IllegalStateException("Cannot add after iterator() has been called");
        }
        ram.add(record);
        if (ram.size() >= maxRecordsInRam) {
            spill();
        }
    }

    /** @return number of temp files written so far. */
    public int getNumSpillFiles() {
        return files.size();
    }

    /** Sorts the chunk in memory and merges equal records in place. */
    private void sortAndCombineRam() {
        ram.sort(comparator);
        int out = 0;
        T current = null;
        for (final T t : ram) {
            if (current != null && comparator.compare(current, t) == 0) {
                current = combiner.apply(current, t);
            } else {
                if (current != null) {
                    ram.set(out++, current);
                }
                current = t;
            }
        }
        if (current != null) {
            ram.set(out++, current);
        }
        ram.subList(out, ram.size()).clear();
    }

    private void spill() {
        sortAndCombineRam();
        try {
            // The htsjdk default temp dir is <java.io.tmpdir>/<user.name>.  Nothing guarantees it exists: it can be
            // gone by the time we spill (e.g. after an earlier htsjdk SortingCollection finished), so create it.
            final File tmpDir = IOUtil.getDefaultTmpDir();
            Files.createDirectories(tmpDir.toPath());
            final File file = IOUtil.newTempFile("CombiningSortingCollection.", ".tmp", new File[]{tmpDir});
            file.deleteOnExit();
            files.add(file);
            try (OutputStream os = tempStreamFactory.wrapTempOutputStream(new FileOutputStream(file), Defaults.BUFFER_SIZE)) {
                codec.setOutputStream(os);
                for (final T t : ram) {
                    codec.encode(t);
                }
            }
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
        ram.clear();
    }

    /** One sorted stream (a spilled chunk or the in-memory remainder) taking part in the merge. */
    private class Source implements Comparable<Source> {
        final int index;
        final Iterator<T> iterator;
        final InputStream stream;
        final File file;
        T head;

        Source(final int index, final Iterator<T> iterator, final InputStream stream, final File file) {
            this.index = index;
            this.iterator = iterator;
            this.stream = stream;
            this.file = file;
            advance();
        }

        /** @return false when exhausted. */
        boolean advance() {
            if (iterator.hasNext()) {
                head = iterator.next();
                return true;
            }
            head = null;
            return false;
        }

        void close() {
            try {
                if (stream != null) {
                    stream.close();
                }
            } catch (IOException e) {
                throw new RuntimeIOException(e);
            } finally {
                if (file != null) {
                    file.delete();
                }
            }
        }

        @Override
        public int compareTo(final Source other) {
            final int cmp = comparator.compare(head, other.head);
            return cmp != 0 ? cmp : Integer.compare(index, other.index);
        }
    }

    private class FileRecordIterator implements Iterator<T> {
        private final SortingCollection.Codec<T> fileCodec = codec.clone();
        private T next;

        FileRecordIterator(final InputStream stream) {
            fileCodec.setInputStream(stream);
            next = fileCodec.decode();
        }

        @Override
        public boolean hasNext() {
            return next != null;
        }

        @Override
        public T next() {
            if (next == null) {
                throw new NoSuchElementException();
            }
            final T ret = next;
            next = fileCodec.decode();
            return ret;
        }
    }

    private class MergingIterator implements CloseableIterator<T> {
        private final PriorityQueue<Source> queue = new PriorityQueue<>();
        private final List<Source> all = new ArrayList<>();

        MergingIterator() {
            try {
                for (int i = 0; i < files.size(); i++) {
                    final File file = files.get(i);
                    final InputStream stream = tempStreamFactory.wrapTempInputStream(new FileInputStream(file), Defaults.BUFFER_SIZE);
                    final Source s = new Source(i, new FileRecordIterator(stream), stream, file);
                    all.add(s);
                    if (s.head != null) {
                        queue.add(s);
                    }
                }
                final Source ramSource = new Source(files.size(), ram.iterator(), null, null);
                all.add(ramSource);
                if (ramSource.head != null) {
                    queue.add(ramSource);
                }
            } catch (IOException e) {
                close();
                throw new RuntimeIOException(e);
            }
        }

        @Override
        public boolean hasNext() {
            return !queue.isEmpty();
        }

        private void advanceAndRequeue(final Source s) {
            if (s.advance()) {
                queue.add(s);
            } else {
                s.close();
            }
        }

        @Override
        public T next() {
            if (queue.isEmpty()) {
                throw new NoSuchElementException();
            }
            final Source first = queue.poll();
            T result = first.head;
            advanceAndRequeue(first);
            while (!queue.isEmpty() && comparator.compare(result, queue.peek().head) == 0) {
                final Source s = queue.poll();
                result = combiner.apply(result, s.head);
                advanceAndRequeue(s);
            }
            return result;
        }

        @Override
        public void close() {
            queue.clear();
            for (final Source s : all) {
                s.close();
            }
            all.clear();
        }
    }

    /**
     * Ends the adding phase. Can be called once.
     * @return the combined records in comparator order. Close it to delete temp files early.
     */
    public CloseableIterator<T> iterator() {
        if (iterating) {
            throw new IllegalStateException("iterator() can only be called once");
        }
        iterating = true;
        sortAndCombineRam();
        if (files.isEmpty()) {
            final Iterator<T> it = ram.iterator();
            return new CloseableIterator<T>() {
                @Override
                public void close() {}

                @Override
                public boolean hasNext() {
                    return it.hasNext();
                }

                @Override
                public T next() {
                    return it.next();
                }
            };
        }
        return new MergingIterator();
    }
}
