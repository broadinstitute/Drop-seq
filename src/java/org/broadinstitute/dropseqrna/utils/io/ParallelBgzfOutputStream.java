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
package org.broadinstitute.dropseqrna.utils.io;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.IOUtil;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayDeque;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.function.IntConsumer;

/**
 * Writes a BGZF file, compressing blocks in parallel.  It is built only on htsjdk's public API: each block of
 * uncompressed data is compressed by its own {@link BlockCompressedOutputStream}, and the file ends with the
 * standard BGZF end-of-file block, which that class writes.  The result is a valid BGZF file (so tools such as
 * tabix and htsjdk understand it), which is also valid gzip: any gzip reader (zcat, java.util.zip.GZIPInputStream,
 * R, Python) reads it as the original bytes.  The compressed bytes differ from those of a single-member
 * GZIPOutputStream, the decompressed data does not.
 * <p>
 * Compression uses htsjdk's process-wide default deflater factory, which the Picard options --USE_JDK_DEFLATER and
 * --USE_JDK_INFLATER set.  With --USE_JDK_DEFLATER true the factory is libdeflate, much faster than JDK zlib; with the
 * Picard default of false it is the Intel deflater if available, otherwise JDK zlib.
 * <p>
 * With {@code numThreads == 1} each block is compressed on the calling thread (no extra threads).  With more
 * threads, blocks are compressed by a pool of {@code numThreads} daemon threads and written in submission order by
 * the calling thread; at most {@code 2 * numThreads} blocks are in flight, which bounds memory and applies
 * backpressure.  A failure in a worker is rethrown from the next {@link #write} or from {@link #close}.
 * <p>
 * Not thread safe: use from one thread at a time, like any OutputStream.  {@link #close()} is required to finish the file.
 */
public class ParallelBgzfOutputStream extends OutputStream {

    public static final int DEFAULT_BLOCK_SIZE = 1024 * 1024;

    private final OutputStream out;
    private final int level;
    private final int blockSize;
    private final int maxInFlight;
    private final ExecutorService pool;                          // null when numThreads == 1
    private final ArrayDeque<Future<byte[]>> inFlight = new ArrayDeque<>();

    private byte[] block;
    private int blockLength = 0;
    private int blocksSubmitted = 0;
    private IOException failure = null;
    private boolean closed = false;

    /** Test hook: called with the block index on the worker just before the block is compressed. */
    volatile IntConsumer beforeBlockDeflate = null;

    /**
     * @param out destination; closed by {@link #close()}
     * @param level compression level, 0 to 9 as for Picard's COMPRESSION_LEVEL
     * @param numThreads number of compression threads; 1 compresses on the calling thread
     */
    public ParallelBgzfOutputStream(final OutputStream out, final int level, final int numThreads) {
        this(out, level, numThreads, DEFAULT_BLOCK_SIZE);
    }

    public ParallelBgzfOutputStream(final OutputStream out, final int level, final int numThreads, final int blockSize) {
        if (numThreads < 1) throw new IllegalArgumentException("numThreads must be >= 1: " + numThreads);
        if (blockSize < 1) throw new IllegalArgumentException("blockSize must be >= 1: " + blockSize);
        this.out = out;
        this.level = level;
        this.blockSize = blockSize;
        this.maxInFlight = 2 * numThreads;
        this.block = new byte[blockSize];
        if (numThreads == 1) {
            pool = null;
        } else {
            pool = Executors.newFixedThreadPool(numThreads, r -> {
                final Thread t = new Thread(r, "ParallelBgzfOutputStream-compress");
                t.setDaemon(true);
                return t;
            });
        }
    }

    /**
     * Opens {@code file} for writing.  A name ending in .gz gets a ParallelBgzfOutputStream at the current Picard
     * compression level; any other name is opened as IOUtil.openFileForWriting does.
     */
    public static OutputStream openFileForWriting(final File file, final int numThreads) {
        if (!file.getName().endsWith(".gz")) {
            return IOUtil.openFileForWriting(file);
        }
        try {
            return new ParallelBgzfOutputStream(Files.newOutputStream(file.toPath()), IOUtil.getCompressionLevel(), numThreads);
        } catch (IOException e) {
            throw new htsjdk.samtools.util.RuntimeIOException("Error opening " + file.getAbsolutePath(), e);
        }
    }

    @Override
    public void write(final int b) throws IOException {
        checkOpen();
        if (blockLength == blockSize) submitBlock();
        block[blockLength++] = (byte) b;
    }

    @Override
    public void write(final byte[] bytes, int offset, int length) throws IOException {
        checkOpen();
        while (length > 0) {
            if (blockLength == blockSize) submitBlock();
            final int n = Math.min(length, blockSize - blockLength);
            System.arraycopy(bytes, offset, block, blockLength, n);
            blockLength += n;
            offset += n;
            length -= n;
        }
    }

    /** Does not force out a partial block: that would make tiny blocks. Data is complete after close(). */
    @Override
    public void flush() throws IOException {
        checkOpen();
    }

    @Override
    public void close() throws IOException {
        if (closed) return;
        closed = true;
        try {
            rethrowFailure();
            if (blockLength > 0) submitBlock();
            while (!inFlight.isEmpty()) writeOldest();
            // The end-of-file block.  With nothing written this is also all that an empty file contains.
            out.write(compress(new byte[0], 0, true));
            out.flush();
        } finally {
            if (pool != null) {
                pool.shutdownNow();
                try {
                    pool.awaitTermination(60, TimeUnit.SECONDS);
                } catch (InterruptedException e) {
                    Thread.currentThread().interrupt();
                }
            }
            out.close();
        }
        rethrowFailure();
    }

    private void checkOpen() throws IOException {
        if (closed) throw new IOException("Stream closed");
        rethrowFailure();
    }

    private void rethrowFailure() throws IOException {
        if (failure != null) throw new IOException("Error compressing gzip block", failure);
    }

    private void submitBlock() throws IOException {
        final byte[] data = block;
        final int length = blockLength;
        final int index = blocksSubmitted++;
        block = new byte[blockSize];
        blockLength = 0;
        if (pool == null) {
            try {
                out.write(compressBlock(data, length, index));
            } catch (IOException e) {
                failure = e;
                throw e;
            }
        } else {
            inFlight.add(pool.submit(() -> compressBlock(data, length, index)));
            while (inFlight.size() > maxInFlight) writeOldest();
        }
    }

    private void writeOldest() throws IOException {
        try {
            out.write(inFlight.poll().get());
        } catch (ExecutionException e) {
            failure = e.getCause() instanceof IOException ? (IOException) e.getCause() : new IOException(e.getCause());
            throw new IOException("Error compressing gzip block", failure);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            failure = new IOException("Interrupted while compressing gzip block", e);
            throw new IOException(failure);
        } catch (IOException e) {
            failure = e;
            throw e;
        }
    }

    private byte[] compressBlock(final byte[] data, final int length, final int index) throws IOException {
        final IntConsumer hook = beforeBlockDeflate;
        if (hook != null) hook.accept(index);
        try {
            return compress(data, length, false);
        } catch (RuntimeException e) {
            throw new IOException("Error compressing gzip block " + index, e);
        }
    }

    /** Compresses one block with htsjdk, with the end-of-file block after it only if {@code endOfFile}. */
    private byte[] compress(final byte[] data, final int length, final boolean endOfFile) throws IOException {
        final ByteArrayOutputStream bytes = new ByteArrayOutputStream();
        // No file: the output goes to the byte array, and the Path is only used by htsjdk to name the stream.
        final BlockCompressedOutputStream bgzf = new BlockCompressedOutputStream(bytes, (Path) null, level);
        try {
            bgzf.write(data, 0, length);
        } finally {
            // Closed exactly once: close() after close(boolean) would write a second end-of-file block.
            bgzf.close(endOfFile);
        }
        return bytes.toByteArray();
    }
}
