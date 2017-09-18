/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2015 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.dropseqrna.utils.readiterators;

import com.google.common.collect.Interner;
import com.google.common.collect.Interners;
import htsjdk.samtools.*;
import htsjdk.samtools.util.IOUtil;
import picard.PicardException;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

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
        final List<SamReader> readers = new ArrayList<SamReader>(inputs.size());
        final List<SAMFileHeader> headers = new ArrayList<SAMFileHeader>(inputs.size());
        final Interner<SAMSequenceDictionary> sequenceDictionaryInterner =Interners.newStrongInterner();
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
                            header.getSortOrder(), inFile.getAbsolutePath(), inputSortOrder, inputs.get(0).getAbsolutePath()));
                }
            }
            headers.add(header);
        }

        if (inputs.size() == 1) {
            return new SamHeaderAndIterator(headers.get(0), readers.get(0).iterator());
        } else {
            final SAMFileHeader.SortOrder outputSortOrder;
            if (maintainSort) {
                outputSortOrder = inputSortOrder;
            } else {
                outputSortOrder = SAMFileHeader.SortOrder.unsorted;
            }
            final SamFileHeaderMerger headerMerger = new SamFileHeaderMerger(outputSortOrder, headers, false);
            final MergingSamRecordIterator iterator = new MergingSamRecordIterator(headerMerger, readers, true);
            return new SamHeaderAndIterator(headerMerger.getMergedHeader(), iterator);
        }
    }

    /**
     *
     * @param maintainSort If true, all inputs must be sorted the same way, and they are merge sorted.  If false, inputs
     *                     are merged in arbitrary order.
     */
    public static SamHeaderAndIterator mergeInputs(final List<File> inputs,
                                                   final boolean maintainSort) {
        return mergeInputs(inputs, maintainSort, SamReaderFactory.makeDefault());
    }
}
