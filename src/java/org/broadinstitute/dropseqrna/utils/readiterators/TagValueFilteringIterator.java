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
 * Filters out a read if the tag value does not match an one of an expected set of values or is not set.
 * If the expected set is null or empty, no reads are filtered.
 * @author nemesh
 *
 */
public class TagValueFilteringIterator<T> extends FilteredIterator<SAMRecord> {
	final short requiredTag;
	final Set<T> expectedValues;
	private static final Log log = Log.getInstance(TagValueFilteringIterator.class);
	
    public TagValueFilteringIterator(final Iterator<SAMRecord> underlyingIterator, final String requiredTag, final Collection<T> expectedValues) {
        super(underlyingIterator);
        this.requiredTag = SAMTag.makeBinaryTag(requiredTag);
        this.expectedValues = new HashSet<T>(expectedValues);
    }



    @Override
    public boolean filterOut(final SAMRecord rec) {
    	Object value = rec.getAttribute(requiredTag);
    	if (value == null)
			return true;

    	if (this.expectedValues.contains(value)) 
    		return false;
        return true;
    }

	@Override
	public void logFilterResults() {		
		String msg = String.format("Required Tag [%s] Records pass [%d] records fail [%d] ",SAMTag.makeStringTag(requiredTag), this.getRecordsPassed(), this.getRecordsFailed());  
		log.info(msg);		
	}
}



