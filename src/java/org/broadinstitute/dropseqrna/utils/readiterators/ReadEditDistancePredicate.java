package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.function.Predicate;

import htsjdk.samtools.SAMRecord;

public class ReadEditDistancePredicate implements Predicate<SAMRecord> {

	private Integer maxEditDistance;
	public static final String EDIT_DISTANCE_TAG="NM";

	public ReadEditDistancePredicate(final Integer maxEditDistance) {
		this.maxEditDistance=maxEditDistance;
	}
	
	@Override
	/**
	 * Returns true when the read has <= threshold edit distance to the reference genome.
	 */
	public boolean test(SAMRecord r) {
		if (maxEditDistance==null) return true;
		Object o = r.getAttribute(EDIT_DISTANCE_TAG);
		if (o==null) return true;
		int readEditDistance = ((Integer) o).intValue();
    	return readEditDistance <= this.maxEditDistance;
	}

}
