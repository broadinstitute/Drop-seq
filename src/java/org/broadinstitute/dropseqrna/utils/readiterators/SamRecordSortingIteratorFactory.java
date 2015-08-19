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

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriterImpl;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;
import org.broadinstitute.dropseqrna.utils.SortingIteratorFactory;

import java.util.Comparator;
import java.util.Iterator;

public class SamRecordSortingIteratorFactory {

    /**
     * @param progressLogger pass null if not interested in progress.
     * @return An iterator with all the records from underlyingIterator, in order defined by comparator.
     */
    public static CloseableIterator<SAMRecord> create(final SAMFileHeader header,
                                           final Iterator<SAMRecord> underlyingIterator,
                                           final Comparator<SAMRecord> comparator,
                                           final ProgressLogger progressLogger) {
        final SortingIteratorFactory.ProgressCallback progressCallback;
        if (progressLogger != null) {
            progressCallback = new SortingIteratorFactory.ProgressCallback<SAMRecord>() {
                @Override
                public void logProgress(SAMRecord record) {
                    progressLogger.record(record);
                }
            };
        } else {
            progressCallback = null;
        }
        return SortingIteratorFactory.create(SAMRecord.class,
                underlyingIterator, comparator, new BAMRecordCodec(header),
                SAMFileWriterImpl.getDefaultMaxRecordsInRam(),
                progressCallback);
    }
}
