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

import htsjdk.samtools.util.BlockCompressedInputStream;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.file.Files;
import java.util.Random;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPInputStream;

/**
 * The output is read back with the JDK's GZIPInputStream, the way any consumer would, and checked to be real BGZF
 * with the standard end-of-file block.
 */
public class ParallelBgzfOutputStreamTest {

    /** Small blocks, so that a few KB of data is many blocks. */
    private static final int BLOCK = 1000;

    /** Compressible text that differs from block to block, so a swapped or lost block changes the result. */
    private static byte[] makeData(final int length) {
        final Random random = new Random(length);
        final StringBuilder sb = new StringBuilder(length + 64);
        while (sb.length() < length) {
            sb.append("CELL").append(random.nextInt(500)).append('\t').append("GENE").append(random.nextInt(50)).append('\n');
        }
        return sb.substring(0, length).getBytes();
    }

    private static byte[] gunzip(final File file) throws IOException {
        try (InputStream in = new GZIPInputStream(new FileInputStream(file))) {
            return in.readAllBytes();
        }
    }

    private static File tempGz() throws IOException {
        final File file = File.createTempFile("ParallelBgzfOutputStreamTest", ".txt.gz");
        file.deleteOnExit();
        return file;
    }

    @DataProvider(name = "numThreads")
    public Object[][] numThreads() {
        return new Object[][]{{1}, {4}};
    }

    /** Many blocks, written as a mix of single bytes and arrays that straddle block boundaries. */
    @Test(dataProvider = "numThreads")
    public void testFileRoundTrip(final int numThreads) throws IOException {
        final byte[] data = makeData(10 * BLOCK + 7);
        final File file = tempGz();
        try (OutputStream out = new ParallelBgzfOutputStream(Files.newOutputStream(file.toPath()), 5, numThreads, BLOCK)) {
            int i = 0;
            while (i < 2500) out.write(data[i++]);
            out.write(data, i, 1234);
            i += 1234;
            out.write(data, i, data.length - i);
        }
        Assert.assertEquals(gunzip(file), data);
        Assert.assertEquals(BlockCompressedInputStream.checkTermination(file), BlockCompressedInputStream.FileTermination.HAS_TERMINATOR_BLOCK);
    }

    /** The public entry point other tools would use: default block size, several blocks, .gz and plain names. */
    @Test
    public void testOpenFileForWriting() throws IOException {
        final byte[] data = makeData(3 * ParallelBgzfOutputStream.DEFAULT_BLOCK_SIZE + 17);
        final File gz = tempGz();
        final File plain = File.createTempFile("ParallelBgzfOutputStreamTest", ".txt");
        plain.deleteOnExit();
        try (OutputStream out = ParallelBgzfOutputStream.openFileForWriting(gz, 3)) {
            out.write(data);
        }
        try (OutputStream out = ParallelBgzfOutputStream.openFileForWriting(plain, 3)) {
            out.write(data);
        }
        Assert.assertEquals(gunzip(gz), data);
        Assert.assertEquals(Files.readAllBytes(plain.toPath()), data, "A name that does not end in .gz must not be compressed");
    }

    /** An empty output file must still be valid gzip (just the end-of-file block). */
    @Test
    public void testEmptyStreamIsValidGzip() throws IOException {
        final File file = tempGz();
        new ParallelBgzfOutputStream(Files.newOutputStream(file.toPath()), 5, 2).close();
        Assert.assertEquals(gunzip(file).length, 0);
        Assert.assertEquals(BlockCompressedInputStream.checkTermination(file), BlockCompressedInputStream.FileTermination.HAS_TERMINATOR_BLOCK);
    }

    /** Block 0 is held until every later block has been compressed, so the blocks finish in reverse order. */
    @Test
    public void testBlocksAreWrittenInOrderWhenTheyFinishOutOfOrder() throws Exception {
        final byte[] data = makeData(4 * BLOCK);
        final File file = tempGz();
        final CountDownLatch laterBlocksDone = new CountDownLatch(3);
        try (ParallelBgzfOutputStream out = new ParallelBgzfOutputStream(Files.newOutputStream(file.toPath()), 5, 4, BLOCK)) {
            out.beforeBlockDeflate = index -> {
                if (index == 0) {
                    try {
                        laterBlocksDone.await(10, TimeUnit.SECONDS);
                    } catch (InterruptedException e) {
                        throw new RuntimeException(e);
                    }
                } else {
                    laterBlocksDone.countDown();
                }
            };
            out.write(data, 0, data.length);
        }
        Assert.assertEquals(laterBlocksDone.getCount(), 0, "Test did not force out of order completion");
        Assert.assertEquals(gunzip(file), data);
    }

    /** A failure on a worker must surface as an IOException (not hang, not silently truncate the file). */
    @Test
    public void testWorkerFailureIsReported() {
        final byte[] data = makeData(3 * BLOCK);
        final ParallelBgzfOutputStream out = new ParallelBgzfOutputStream(new ByteArrayOutputStream(), 5, 2, BLOCK);
        out.beforeBlockDeflate = index -> {
            if (index == 1) throw new IllegalStateException("boom");
        };
        Assert.assertThrows(IOException.class, () -> {
            out.write(data, 0, data.length);
            out.close();
        });
    }
}
