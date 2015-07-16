package org.broadinstitute.dropseqrna.utils.bamtagcomparator;

import htsjdk.samtools.util.Interval;

public class IntervalComparator implements TagComparator {

	public IntervalComparator () {
		
	}
	
	@Override
	public int compare(Object valueOne, Object valueTwo) {
		// TODO Auto-generated method stub
		return 0;
	}
	
	
	
	
	Interval parseInterval (Object val) {
		if (val == null) return null;
		String valS=null;
		if (val instanceof String) {
			valS=(String) val;
		}
		
		return null;
	}

}
