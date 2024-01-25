package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.dropseqrna.utils.PredicateFilteredIterator;

import java.util.Iterator;

public class IntervalFilteringIterator extends PredicateFilteredIterator<SAMRecord> {

    public IntervalFilteringIterator (final Iterator<SAMRecord> underlyingIterator, final IntervalList intervals, final boolean retainIntervals) {
        super(underlyingIterator, new IntervalListPredicate(intervals, retainIntervals));
    }

}
