/*
 * MIT License
 *
 * Copyright 2017 Broad Institute
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
package org.broadinstitute.dropseqrna.utils.readiterators;

import com.google.common.collect.Interner;
import com.google.common.collect.Interners;
import htsjdk.samtools.*;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.dropseqrna.utils.CustomBAMIterators;
import org.broadinstitute.dropseqrna.utils.FileUtils;
import org.broadinstitute.dropseqrna.utils.readiterators.parallel.BamDecodePool;
import org.broadinstitute.dropseqrna.utils.readiterators.parallel.ReadMode;
import picard.PicardException;

import java.io.File;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.BiFunction;
import java.util.function.Function;

/**
 * Utilities for creating a single stream of SAMRecords from multiple input files.
 */
public class SamFileMergeUtil {

    /**
     * Use this overload if you want to override defaults in SamReaderFactory, e.g. eager decode
     * @param maintainSort If true, all inputs must be sorted the same way, and they are merge sorted.  If false, inputs
     *                     are merged in arbitrary order.
     */
    public static SamHeaderAndIterator mergeInputs(final List<File> inputs,
                                                   final boolean maintainSort,
                                                   final SamReaderFactory samReaderFactory) {
        return mergeInputPaths(IOUtil.filesToPaths(inputs), maintainSort, samReaderFactory);
    }

    /**
     * Use this overload if you want to override defaults in SamReaderFactory, e.g. eager decode
     * @param maintainSort If true, all inputs must be sorted the same way, and they are merge sorted.  If false, inputs
     *                     are merged in arbitrary order.
     */
    public static SamHeaderAndIterator mergeInputPaths(final List<Path> inputs,
                                                       final boolean maintainSort,
                                                       final SamReaderFactory samReaderFactory) {
        return mergeInputPaths(inputs, maintainSort, samReaderFactory::open);
    }

    /**
     * Use this overload to read and decode BAM records on multiple threads.  Same as
     * {@link #mergeInputPaths(List, boolean, SamReaderFactory, int, ReadMode)} with {@link ReadMode#AUTO}.
     * @param maintainSort If true, all inputs must be sorted the same way, and they are merge sorted.  If false, inputs
     *                     are merged in arbitrary order.
     * @param numThreads Number of threads that read and decode BAM records, shared by all inputs.  1 is the same as
     *                   the overload without this argument.
     *                   <p>
     *                   Whether more threads help depends on the input (from benchmarks; check it against your own
     *                   data):
     *                   <ul>
     *                   <li>Several inputs: reading scales with many threads.  They are read at the same time, so
     *                   the speedup follows the number of inputs, roughly one thread per input.</li>
     *                   <li>One aligned input: reading scales with up to 2 threads.  An uncompressed input kept
     *                   speeding up to about 4 threads.</li>
     *                   <li>One unaligned input: asynchronous IO is better.  Such reads are cheap to decode, so
     *                   extra threads help little.  With more than 1 thread this input is read with htsjdk's
     *                   asynchronous IO instead, see {@link ReadMode#AUTO}, and numThreads does not set its number
     *                   of threads.</li>
     *                   </ul>
     *                   <p>
     *                   With more than 1 thread, several inputs and maintainSort false, all inputs are read at the
     *                   same time and records come one batch from each input in turn, not in file order.  The order
     *                   depends on numThreads being above 1 but is otherwise deterministic.  Callers that depend on
     *                   record order must not combine this with maintainSort false.
     */
    public static SamHeaderAndIterator mergeInputPaths(final List<Path> inputs,
                                                       final boolean maintainSort,
                                                       final SamReaderFactory samReaderFactory,
                                                       final int numThreads) {
        return mergeInputPaths(inputs, maintainSort, samReaderFactory, numThreads, ReadMode.AUTO);
    }

    /**
     * As {@link #mergeInputPaths(List, boolean, SamReaderFactory, int)}, with a choice of how the inputs are read.
     * See {@link ReadMode}, and {@link ReadMode#resolve} for what is used.  With {@link ReadMode#STOCK} or
     * {@link ReadMode#ASYNC}, and several unsorted inputs, the inputs are read one after another in file order.
     */
    public static SamHeaderAndIterator mergeInputPaths(final List<Path> inputs,
                                                       final boolean maintainSort,
                                                       final SamReaderFactory samReaderFactory,
                                                       final int numThreads,
                                                       final ReadMode mode) {
        switch (mode.resolve(inputs, samReaderFactory, numThreads)) {
            case ASYNC:
                return mergeInputPaths(inputs, maintainSort, asyncFactory(samReaderFactory));
            case STOCK:
                return mergeInputPaths(inputs, maintainSort, samReaderFactory);
            default:
                final BamDecodePool pool = new BamDecodePool(numThreads, samReaderFactory);
                final SamHeaderAndIterator merged = mergeInputPaths(inputs, maintainSort, pool::open, pool::interleave);
                return new SamHeaderAndIterator(merged.header, pool.closeWhenClosed(merged.iterator));
        }
    }

    /** A new factory, so that the caller's is not changed: eager decode, async IO, the caller's validation. */
    private static SamReaderFactory asyncFactory(final SamReaderFactory samReaderFactory) {
        return SamReaderFactory.makeDefault()
                .validationStringency(samReaderFactory.validationStringency())
                .enable(SamReaderFactory.Option.EAGERLY_DECODE)
                .setUseAsyncIo(true);
    }

    /**
     * Use this overload to control how each input is opened, e.g. with a thread pool that decodes records in parallel.
     * @param maintainSort If true, all inputs must be sorted the same way, and they are merge sorted.  If false, inputs
     *                     are merged in arbitrary order.
     * @param opener Opens one input.  The readers it returns are closed when the returned iterator is closed.
     */
    public static SamHeaderAndIterator mergeInputPaths(final List<Path> inputs,
                                                       final boolean maintainSort,
                                                       final Function<Path, SamReader> opener) {
        return mergeInputPaths(inputs, maintainSort, opener, UnsortedMergingSamRecordIterator::new);
    }

    private static SamHeaderAndIterator mergeInputPaths(final List<Path> inputs,
                                                        final boolean maintainSort,
                                                        final Function<Path, SamReader> opener,
                                                        final BiFunction<SAMFileHeader, List<SamReader>, CloseableIterator<SAMRecord>> unsortedMerger) {
        if (inputs.isEmpty()) {
            throw new IllegalArgumentException("At least one input must be provided");
        }
        final List<SamReader> readers = new ArrayList<>(inputs.size());
        final List<SAMFileHeader> headers = new ArrayList<>(inputs.size());
        final Interner<SAMSequenceDictionary> sequenceDictionaryInterner = Interners.newStrongInterner();
        SAMFileHeader.SortOrder inputSortOrder = null;
        for (final Path inFile : inputs) {
            IOUtil.assertFileIsReadable(inFile);
            final SamReader in = opener.apply(inFile);
            readers.add(in);
            final SAMFileHeader header = in.getFileHeader();
            header.setSequenceDictionary(sequenceDictionaryInterner.intern(header.getSequenceDictionary()));
            if (maintainSort) {
                if (inputSortOrder == null) {
                    inputSortOrder = header.getSortOrder();
                } else if (header.getSortOrder() != inputSortOrder) {
                    throw new PicardException(String.format("Sort order(%s) of %s does not agree with sort order(%s) of %s",
                            header.getSortOrder(), FileUtils.toAbsoluteString(inFile), inputSortOrder,
                            FileUtils.toAbsoluteString(inputs.getFirst())));
                }
            }
            headers.add(header);
        }

        if (inputs.size() == 1) {
            return new SamHeaderAndIterator(headers.getFirst(), readers.getFirst().iterator());
        } else {
            final SAMFileHeader.SortOrder outputSortOrder;
            if (maintainSort) {
                outputSortOrder = inputSortOrder;
            } else {
                outputSortOrder = SAMFileHeader.SortOrder.unsorted;
            }
            final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(outputSortOrder, headers, false);
            final CloseableIterator<SAMRecord> iterator;
            if (outputSortOrder == SAMFileHeader.SortOrder.unsorted) {
                // If unsorted, output the reads in the same order they came in, unless the caller supplied a
                // different merger.
                // In contrast, MergingSamRecordIterator re-orders in an arbitrary way based on
                // identityHashCode.
                iterator = unsortedMerger.apply(headerMerger.getMergedHeader(), readers);
            } else {
                iterator = new MergingSamRecordIterator(headerMerger, readers, true);
            }
            return new SamHeaderAndIterator(headerMerger.getMergedHeader(), iterator);
        }
    }

    /**
     * @param maintainSort If true, all inputs must be sorted the same way, and they are merge sorted.  If false, inputs
     *                     are merged in arbitrary order.
     */
    public static SamHeaderAndIterator mergeInputs(final List<File> inputs,
                                                   final boolean maintainSort) {
        return mergeInputPaths(IOUtil.filesToPaths(inputs), maintainSort, SamReaderFactory.makeDefault());
    }

    /**
     * @param maintainSort If true, all inputs must be sorted the same way, and they are merge sorted.  If false, inputs
     *                     are merged in arbitrary order.
     */
    public static SamHeaderAndIterator mergeInputPaths(final List<Path> inputs,
                                                       final boolean maintainSort) {
        return mergeInputPaths(inputs, maintainSort, SamReaderFactory.makeDefault());
    }

    /**
     * Read a collection of BAMs in the desired order.  Headers must have compatible sequence dictionaries.
     * @param inputs one or more BAMs, in any order.
     * @param sortOrder Desired sort order
     * @param samReaderFactory
     * @return
     */
    public static SamHeaderAndIterator mergeInputs(final List<File> inputs,
                                                   final SAMFileHeader.SortOrder sortOrder,
                                                   final SamReaderFactory samReaderFactory,
                                                   final boolean mergeDictionaries) {
        return mergeInputPaths(IOUtil.filesToPaths(inputs), sortOrder, samReaderFactory, mergeDictionaries);
    }

    /**
     * Read a collection of BAMs in the desired order.  Headers must have compatible sequence dictionaries.
     *
     * @param inputs           one or more BAMs, in any order.
     * @param sortOrder        Desired sort order
     * @param samReaderFactory
     * @return
     */
    public static SamHeaderAndIterator mergeInputPaths(final List<Path> inputs,
                                                       final SAMFileHeader.SortOrder sortOrder,
                                                       final SamReaderFactory samReaderFactory,
                                                       final boolean mergeDictionaries) {
        final List<SamReader> readers = inputs
                .stream()
                .map(input -> samReaderFactory.open(input))
                .toList();
        final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(sortOrder, readers.stream().map(SamReader::getFileHeader).toList(), mergeDictionaries);
        final CloseableIterator<SAMRecord> iterator;
        if (sortOrder == SAMFileHeader.SortOrder.unsorted) {
            iterator = new UnsortedMergingSamRecordIterator(headerMerger.getMergedHeader(), readers);
        } else {
            final Map<SamReader, CloseableIterator<SAMRecord>> map = new HashMap<>(readers.size());
            for (final SamReader reader : readers) {
                map.put(reader,
                        CustomBAMIterators.getSortedRecords(new SamHeaderAndIterator(reader.getFileHeader(), reader.iterator()), sortOrder));
            }
            iterator = new MergingSamRecordIterator(headerMerger, map, true);
        }
        return new SamHeaderAndIterator(headerMerger.getMergedHeader(), iterator);
    }
}
