package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;

import java.util.Iterator;

import org.broadinstitute.dropseqrna.utils.FilteredIterator;

public class PCRDuplicateFilteringIterator extends FilteredIterator<SAMRecord> {

	public PCRDuplicateFilteringIterator(final Iterator<SAMRecord> underlyingIterator) {
		super(underlyingIterator);

	}

	@Override
	public boolean filterOut(final SAMRecord rec) {
		return rec.getDuplicateReadFlag();

	}

}
