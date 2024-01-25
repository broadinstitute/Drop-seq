package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.OverlapDetector;

import java.util.function.Predicate;

public class IntervalListPredicate implements Predicate<SAMRecord> {
    private final OverlapDetector<Interval> od;
    private final boolean retainIntervals;
    public IntervalListPredicate (IntervalList intervalList, boolean retainIntervals) {
        od = new OverlapDetector<>(0, 0);
        od.addAll(intervalList.getIntervals(), intervalList.getIntervals());
        this.retainIntervals=retainIntervals;
    }

    @Override
    /**
     * Returns true when the read overlaps an interval in the interval list
     */
    public boolean test(SAMRecord rec) {
        Interval i = new Interval (rec.getContig(), rec.getStart(), rec.getEnd());
        boolean overlapFound = od.overlapsAny(i);
        // if the overlap is found and you're retaining, don't filter.  If overlap not found and not retaining, don't filter.
        return overlapFound == retainIntervals;
    }
}
