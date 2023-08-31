package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.Iterator;

import org.broadinstitute.dropseqrna.utils.FilteredIterator;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;

public class OverlapFilteringIterator<OD_TYPE> extends FilteredIterator<SAMRecord> {
	private static final Log log = Log.getInstance(OverlapFilteringIterator.class);
	@SuppressWarnings("rawtypes")
	private final OverlapDetector<OD_TYPE> od;
	private final boolean retainOverlap;
	
	/**
	 * Filter or retain records that overlap intervals in an overlap detector.  
	 * @param underlyingIterator
	 * @param od The set of intervals to check
	 * @param retainOverlap If true, then retain reads that overlap an interval and filter all others.  If false, filter reads that overlap intervals.
	 */
	public OverlapFilteringIterator(final Iterator<SAMRecord> underlyingIterator, final OverlapDetector<OD_TYPE> od, boolean retainOverlap) {
		super(underlyingIterator);
		this.od=od;
		this.retainOverlap=retainOverlap;
	}

	@Override
	/**
	 * For a given record test to see if it overlaps intervals in the overlap detector.
	 * 
	 * @param rec
	 * @return
	 */
	public boolean filterOut(final SAMRecord rec) {		
		if (od.overlapsAny(rec)) {
			// if this read overlaps an interval and retain record is true, then return false. 
			return !retainOverlap;
		}
		// if this read does not overlap and interval and retain record is true, return true to filter the read.
		return retainOverlap;		 
	}

	@Override
	public void logFilterResults() {
		String msg = String.format("Records pass [%d] records fail [%d] ",this.getRecordsPassed(), this.getRecordsFailed());  
		log.info(msg);										
	}


}
