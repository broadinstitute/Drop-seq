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
package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.SortingCollection;

import java.util.Comparator;
import java.util.Iterator;

public class SortingIteratorFactory {

    /**
     * Create one of these if you want to get progress as records are being stuffed into the SortingCollection.
     *
     */
    public interface ProgressCallback<T> {
        public void logProgress(T record);
    }

    /**
     *
     * @param componentType Required because of Java generic syntax limitations.
     * @param underlyingIterator All records are pulled from this iterator, which is then closed if closeable.
     * @param comparator Defines sort order.
     * @param codec For spilling to temp files
     * @param maxRecordsInRam
     * @param progressLogger Pass null if not interested in logging.
     * @return An iterator in the order defined by comparator, that will produce all the records from underlyingIterator.
     */
    public static <T> CloseableIterator<T> create(final Class<T> componentType,
                                                  final Iterator<T> underlyingIterator,
                                                  final Comparator<T> comparator,
                                                  final SortingCollection.Codec<T> codec,
                                                  final int maxRecordsInRam,
                                                  final ProgressCallback progressLogger) {

        SortingCollection<T> sortingCollection =
                SortingCollection.newInstance(componentType, codec, comparator, maxRecordsInRam);

        while (underlyingIterator.hasNext()) {
            final T rec = underlyingIterator.next();
            if (progressLogger != null)
				progressLogger.logProgress(rec);
            sortingCollection.add(rec);
        }
        CloseableIterator<T> ret = sortingCollection.iterator();
        CloserUtil.close(underlyingIterator);
        return ret;
    }
}
