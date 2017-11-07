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

public class MapQualityFilteredIterator extends FilteredIterator<SAMRecord> {

	private final Integer mapQuality;
	private final boolean rejectNonPrimaryReads;

	public MapQualityFilteredIterator(final Iterator<SAMRecord> underlyingIterator, final Integer mapQuality, final boolean rejectNonPrimaryReads) {
        super(underlyingIterator);
		this.mapQuality=mapQuality;
		this.rejectNonPrimaryReads=rejectNonPrimaryReads;
	}

    @Override
    public boolean filterOut(final SAMRecord r) {
    	/*
    	if (r.getReadName().equals("HN7TNBGXX:4:11609:2753:3252"))
			System.out.println("STOP");
		*/
    	boolean flag=(rejectNonPrimaryReads && r.isSecondaryOrSupplementary()) ||
        (this.mapQuality!=null && r.getMappingQuality() < this.mapQuality);
        return flag;
    }
}
