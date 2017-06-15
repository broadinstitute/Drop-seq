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

import java.util.Comparator;
import java.util.List;

/**
 * Create a comparator from a list of comparators that all operate on the same type.  The first comparator is the
 * primary sort, second comparator is secondary sort, etc.
 */
public class MultiComparator<T> implements Comparator<T> {

    private final Comparator<T>[] comparators;

    @SuppressWarnings("unchecked")
    public MultiComparator(final List<Comparator<T>> comparators) {
        this.comparators = comparators.toArray(new Comparator[comparators.size()]);
    }

    public MultiComparator(final Comparator<T>...comparators) {
        this.comparators = comparators;
    }

    @Override
    public int compare(final T o1, final T o2) {
        for (final Comparator<T> comparator : comparators) {
            final int cmp = comparator.compare(o1, o2);
            if (cmp != 0)
				return cmp;
        }
        return 0;
    }
}
