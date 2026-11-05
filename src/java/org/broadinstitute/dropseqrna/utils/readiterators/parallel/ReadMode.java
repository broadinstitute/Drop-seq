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

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;

import java.nio.file.Path;
import java.util.List;

/**
 * How BAM files are read when a number of threads is requested, see
 * {@link org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil}.
 */
public enum ReadMode {
    /**
     * The default.  With numThreads 1 this is {@link #STOCK}.  Otherwise it is {@link #POOL}, except for a single
     * input whose header has an empty sequence dictionary (unaligned reads), which is read as {@link #ASYNC}.  Such
     * an input is cheap to decode, so the pool helps little, and asynchronous IO was about as fast as the pool on
     * uncompressed files and faster on compressed ones.  An input with a sequence dictionary but only unmapped
     * reads counts as aligned.
     */
    AUTO,
    /** The stock htsjdk reader on the calling thread, whatever numThreads is. */
    STOCK,
    /** Read and decode on a {@link BamDecodePool} when numThreads is above 1, whatever the inputs; else {@link #STOCK}. */
    POOL,
    /**
     * htsjdk's asynchronous IO, whatever numThreads is, and never a pool.  It uses its own reader threads, so the
     * thread count is not numThreads.  Records are decoded eagerly, and validated with the validation stringency of
     * the given factory.  Its other options are not applied, as for {@link BamDecodePool}.  The given factory is
     * not changed.
     */
    ASYNC;

    /**
     * @return what will be used to read these inputs with this mode: {@link #STOCK}, {@link #POOL} or {@link #ASYNC},
     * never {@link #AUTO}.  Only {@link #AUTO} with several threads and one input opens a file, to read its header.
     */
    public ReadMode resolve(final List<Path> inputs, final SamReaderFactory samReaderFactory, final int numThreads) {
        switch (this) {
            case STOCK:
            case ASYNC:
                return this;
            case POOL:
                return numThreads > 1 ? POOL : STOCK;
            default:
                if (numThreads == 1) {
                    return STOCK;
                }
                return inputs.size() == 1 && isUnaligned(inputs.getFirst(), samReaderFactory) ? ASYNC : POOL;
        }
    }

    /** @return true if the header of the file has no sequence dictionary, which is what unaligned reads have. */
    private static boolean isUnaligned(final Path input, final SamReaderFactory samReaderFactory) {
        final SamReader reader = samReaderFactory.open(input);
        final boolean unaligned = reader.getFileHeader().getSequenceDictionary().isEmpty();
        CloserUtil.close(reader);
        return unaligned;
    }
}
