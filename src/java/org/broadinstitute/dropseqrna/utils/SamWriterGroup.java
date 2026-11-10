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

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.dropseqrna.utils.batchwriter.BatchSink;
import org.broadinstitute.dropseqrna.utils.batchwriter.BatchSinkGroup;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * A numbered group of output SAM or BAM files that share one header.  Add outputs with {@link #addOutput}, give records
 * to them with {@link #write}, then {@link #close()}.  Create it with one of two factory methods:
 * <ul>
 * <li>{@link #direct}: the calling thread encodes and writes every record itself, in {@link #write}, straight to the
 * htsjdk writer.  No other threads or batching are involved, and a record can be reused as soon as {@link #write}
 * returns.</li>
 * <li>{@link #batched}: {@link #write} only collects the record, and other threads, shared by all the outputs, encode
 * and write it later in batches.  A record must not be changed or reused after it is written.  A failure in a writer
 * thread is thrown by a later {@link #write} or by {@link #close()}.</li>
 * </ul>
 * Either way, only one thread calls the methods of the group, and each output gets its records in the order written.
 * <p>
 * With a single output, {@link #batched} gains little from more than one writer thread: an output's batches are
 * written in order, one at a time, so at most one writer thread is ever busy for it, and the extra ones sit idle.  The
 * calling thread overlapping with that one writer thread is the whole benefit.  Measured on 40M records, 1 output,
 * 2 input threads: 26.4 s with 2 threads and 27.7 s with 4, against 40.4 s with 1 (see
 * {@code ~/split_bam_by_cell_tuning_plan.md}).  It is not slower than htsjdk's async IO there (27.6 s and 30.0 s), so
 * it is not an anti-pattern, but more threads only cost CPU.  Use at most 2 threads for a single output.
 */
public abstract class SamWriterGroup {

    private final SAMFileHeader header;
    private final SAMFileWriterFactory samWriterFactory;
    private final List<Output> outputs = new ArrayList<>();

    /** One output file, and how many records have been written to it. */
    private abstract static class Output {
        long readCount = 0;

        abstract void add(SAMRecord record);
    }

    /** The calling thread writes every record itself. */
    public static SamWriterGroup direct(final SAMFileHeader header, final boolean createIndex) {
        return new Direct(header, createIndex);
    }

    /**
     * @param writerThreads number of threads that encode and write the outputs, shared by all of them, at least 1
     * @param batchSize     records handed to the writer threads at a time
     */
    public static SamWriterGroup batched(final SAMFileHeader header, final boolean createIndex, final int writerThreads,
                                         final int batchSize) {
        return new Batched(header, createIndex, writerThreads, batchSize);
    }

    private SamWriterGroup(final SAMFileHeader header, final boolean createIndex) {
        this.header = header;
        this.samWriterFactory = new SAMFileWriterFactory().setCreateIndex(createIndex);
    }

    /** Creates the file and returns the index of the new output, which is the number of outputs there were before. */
    public int addOutput(final File file) {
        outputs.add(newOutput(samWriterFactory.makeSAMOrBAMWriter(header, true, file)));
        return outputs.size() - 1;
    }

    public int size() {
        return outputs.size();
    }

    public void write(final int outputIdx, final SAMRecord record) {
        final Output output = outputs.get(outputIdx);
        output.add(record);
        output.readCount++;
    }

    /** The number of records written to the output so far. */
    public long getReadCount(final int outputIdx) {
        return outputs.get(outputIdx).readCount;
    }

    /** Writes any records still waiting, and closes all the outputs. */
    public abstract void close();

    abstract Output newOutput(SAMFileWriter writer);

    private static class Direct extends SamWriterGroup {
        private final List<SAMFileWriter> writers = new ArrayList<>();

        Direct(final SAMFileHeader header, final boolean createIndex) {
            super(header, createIndex);
        }

        @Override
        Output newOutput(final SAMFileWriter writer) {
            writers.add(writer);
            return new Output() {
                @Override
                void add(final SAMRecord record) {
                    writer.addAlignment(record);
                }
            };
        }

        @Override
        public void close() {
            writers.forEach(SAMFileWriter::close);
        }
    }

    private static class Batched extends SamWriterGroup {
        private final BatchSinkGroup group;

        Batched(final SAMFileHeader header, final boolean createIndex, final int writerThreads, final int batchSize) {
            super(header, createIndex);
            this.group = new BatchSinkGroup(writerThreads, batchSize, BatchSinkGroup.DEFAULT_MAX_BATCHES_IN_FLIGHT);
        }

        @Override
        Output newOutput(final SAMFileWriter writer) {
            // The group's threads call the writer one batch at a time, in order, so it needs no synchronization.
            final BatchSink<SAMRecord> sink = group.open(batch -> batch.forEach(writer::addAlignment), writer::close);
            return new Output() {
                @Override
                void add(final SAMRecord record) {
                    sink.add(record);
                }
            };
        }

        @Override
        public void close() {
            group.close();
        }
    }
}
