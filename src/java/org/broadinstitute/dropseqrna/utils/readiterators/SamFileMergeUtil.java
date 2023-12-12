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
import picard.PicardException;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
        if (inputs.isEmpty()) {
            throw new IllegalArgumentException("At least one input must be provided");
        }
        final List<SamReader> readers = new ArrayList<>(inputs.size());
        final List<SAMFileHeader> headers = new ArrayList<>(inputs.size());
        final Interner<SAMSequenceDictionary> sequenceDictionaryInterner = Interners.newStrongInterner();
        SAMFileHeader.SortOrder inputSortOrder = null;
        for (final File inFile : inputs) {
            IOUtil.assertFileIsReadable(inFile);
            final SamReader in = samReaderFactory.open(inFile);
            readers.add(in);
            final SAMFileHeader header = in.getFileHeader();
            header.setSequenceDictionary(sequenceDictionaryInterner.intern(header.getSequenceDictionary()));
            if (maintainSort) {
                if (inputSortOrder == null) {
                    inputSortOrder = header.getSortOrder();
                } else if (header.getSortOrder() != inputSortOrder) {
                    throw new PicardException(String.format("Sort order(%s) of %s does not agree with sort order(%s) of %s",
                            header.getSortOrder(), inFile.getAbsolutePath(), inputSortOrder, inputs.getFirst().getAbsolutePath()));
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
                // If unsorted, output the reads in the same order they came in.
                // In contrast, MergingSamRecordIterator re-orders in an arbitrary way based on
                // identityHashCode.
                iterator = new UnsortedMergingSamRecordIterator(headerMerger.getMergedHeader(), readers);
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
        return mergeInputs(inputs, maintainSort, SamReaderFactory.makeDefault());
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
        final List<SamReader> readers = inputs.stream().map(samReaderFactory::open).toList();
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