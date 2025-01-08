package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.broadinstitute.dropseqrna.utils.FilteredIterator;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.util.Log;

/**
 * If retainExpectedValues is true, then filters out a read if the tag value does not match an one of an expected set of values or is not set.
 * If the expected set is null or empty, no reads are filtered.
 *
 * If retainExpectedValues is false, reject reads that are in the expected values.
 *
 * @author nemesh
 *
 */
public class TagValueFilteringIterator<T> extends FilteredIterator<SAMRecord> {
	final short requiredTag;
	final Set<T> expectedValues;
	private final boolean retainExpectedValues;
	private static final Log log = Log.getInstance(TagValueFilteringIterator.class);

	/**
	 * The default behavior is to retain reads where the tag is not set or has a value in the accepted set.
	 * @param underlyingIterator The iterator to filter
	 * @param requiredTag The tag to filter on
	 * @param expectedValues The values to filter on
	 */
	public TagValueFilteringIterator(final Iterator<SAMRecord> underlyingIterator, final String requiredTag, final Collection<T> expectedValues) {
        this(underlyingIterator, requiredTag, expectedValues, true);
    }

	public TagValueFilteringIterator(final Iterator<SAMRecord> underlyingIterator, final String requiredTag, final Collection<T> expectedValues, boolean retainExpectedValues) {
		super(underlyingIterator);
		this.requiredTag = SAMTag.makeBinaryTag(requiredTag);
		this.expectedValues = new HashSet<T>(expectedValues);
		this.retainExpectedValues=retainExpectedValues;
	}

    @Override
    public boolean filterOut(final SAMRecord rec) {
    	Object value = rec.getAttribute(requiredTag);
    	if (value == null)
			return retainExpectedValues;

    	if (this.expectedValues.contains(value)) 
    		return !retainExpectedValues;
        return retainExpectedValues;
    }

	@Override
	public void logFilterResults() {
		if (this.retainExpectedValues) {
			String msg = String.format("Required Tag [%s] Records pass [%d] records fail [%d] ",SAMTag.makeStringTag(requiredTag), this.getRecordsPassed(), this.getRecordsFailed());
			log.info(msg);
			return;
		}
		String msg = String.format("Rejected Tag [%s] Records pass [%d] records fail [%d] ",SAMTag.makeStringTag(requiredTag), this.getRecordsPassed(), this.getRecordsFailed());
		log.info(msg);		
	}
}



