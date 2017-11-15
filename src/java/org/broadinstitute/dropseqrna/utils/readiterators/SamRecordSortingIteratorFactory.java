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

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriterImpl;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.ProgressLogger;

import java.util.Comparator;
import java.util.Iterator;

import org.broadinstitute.dropseqrna.utils.SortingIteratorFactory;

public class SamRecordSortingIteratorFactory {

    /**
     * @param progressLogger pass null if not interested in progress.
     * @return An iterator with all the records from underlyingIterator, in order defined by comparator.
     */
    public static CloseableIterator<SAMRecord> create(final SAMFileHeader header,
                                           final Iterator<SAMRecord> underlyingIterator,
                                           final Comparator<SAMRecord> comparator,
                                           final ProgressLogger progressLogger) {
        final SortingIteratorFactory.ProgressCallback<SAMRecord> progressCallback;
        if (progressLogger != null)
			progressCallback = new SortingIteratorFactory.ProgressCallback<SAMRecord>() {
                @Override
                public void logProgress(final SAMRecord record) {
                    progressLogger.record(record);
                }
            };
		else
			progressCallback = null;
        return SortingIteratorFactory.create(SAMRecord.class,
                underlyingIterator, comparator, new BAMRecordCodec(header),
                SAMFileWriterImpl.getDefaultMaxRecordsInRam(),
                progressCallback);
    }
}
