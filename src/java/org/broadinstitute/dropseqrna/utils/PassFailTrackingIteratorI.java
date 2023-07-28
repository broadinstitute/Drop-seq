package org.broadinstitute.dropseqrna.utils;

/**
 * Interface to tag iterators that track how many reads pass or fail
 * @author nemesh
 *
 */
public interface PassFailTrackingIteratorI {

	public int getRecordsPassed();
	public int getRecordsFailed();
	
}
