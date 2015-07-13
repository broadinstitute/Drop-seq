package org.broadinstitute.dropseqrna.utils.bamtagcomparator;

import java.util.Comparator;

public interface TagComparator extends Comparator {

	public int compare(final Object valueOne, final Object valueTwo);
	
	
}
