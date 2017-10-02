package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.Iterator;

import org.broadinstitute.dropseqrna.utils.FilteredIterator;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

import htsjdk.samtools.SAMRecord;

public class BamTagCountingIterator extends FilteredIterator<SAMRecord> {

	private String tag;
    private ObjectCounter<String> counter;

	public BamTagCountingIterator(final Iterator<SAMRecord> underlyingIterator, final String tag) {
		super(underlyingIterator);
		if (tag!=null) this.counter = new ObjectCounter<>();
		this.tag=tag;
	}

	@Override
	public boolean filterOut(final SAMRecord rec) {
		if (this.tag==null) return false;
		String value = rec.getStringAttribute(tag);
		if (value!=null) counter.increment(value);
		return false;
	}

	public ObjectCounter<String> getCounts () {
		return this.counter;
	}


}
