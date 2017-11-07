/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
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
