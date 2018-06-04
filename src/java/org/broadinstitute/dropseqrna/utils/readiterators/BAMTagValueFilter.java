package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;

import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

public class BAMTagValueFilter extends FilteredIterator<SAMRecord>{

	private String tagName;
	private Set<String> tagValues;

	/**
	 * Filter out reads that have a tag without a value contained in the tagValues collection.
	 * If the tagValues collection is empty or null, then this filter does not filter any records.
	 * @param underlyingIterator
	 * @param tagName
	 * @param tagValues
	 */
	protected BAMTagValueFilter(final Iterator<SAMRecord> underlyingIterator, final String tagName, final Collection<String> tagValues) {
		super(underlyingIterator);
		this.tagName=tagName;
		this.tagValues=new HashSet<String>(tagValues);
	}

	@Override
	public boolean filterOut(final SAMRecord rec) {
		if (tagValues==null) return false;
		if (tagValues.isEmpty()) return false;
		if (tagValues.contains(rec.getStringAttribute(this.tagName))) return false;
		return true;
	}



}
