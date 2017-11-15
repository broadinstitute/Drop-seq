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
