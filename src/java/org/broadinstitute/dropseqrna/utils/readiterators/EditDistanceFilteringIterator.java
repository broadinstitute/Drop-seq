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

import java.util.Iterator;

import org.broadinstitute.dropseqrna.utils.FilteredIterator;

import htsjdk.samtools.SAMRecord;

public class EditDistanceFilteringIterator extends FilteredIterator<SAMRecord> {

	private Integer maxEditDistance;
	private static final String EDIT_DISTANCE_TAG="NM";

	public EditDistanceFilteringIterator (final Iterator<SAMRecord> underlyingIterator, final Integer maxEditDistance) {
		super(underlyingIterator);
		this.maxEditDistance=maxEditDistance;
	}

	@Override
	public boolean filterOut(final SAMRecord r) {
		if (maxEditDistance==null) return false;
		Object o = r.getAttribute(EDIT_DISTANCE_TAG);
		if (o==null) return false;
		int readEditDistance = ((Integer) o).intValue();
    	return readEditDistance > this.maxEditDistance;
    }
}
