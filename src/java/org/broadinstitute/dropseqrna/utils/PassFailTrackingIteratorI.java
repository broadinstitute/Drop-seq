package org.broadinstitute.dropseqrna.utils;

/**
 * Interface to tag iterators that track how many reads pass or fail
 * @author nemesh
 *
 */
public interface PassFailTrackingIteratorI {

	public long getRecordsPassed();
	public long getRecordsFailed();
	
}
