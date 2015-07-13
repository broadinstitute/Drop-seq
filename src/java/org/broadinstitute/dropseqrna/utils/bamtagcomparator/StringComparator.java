package org.broadinstitute.dropseqrna.utils.bamtagcomparator;

import htsjdk.samtools.SAMException;

public class StringComparator implements TagComparator {

	private Short tagName;
	
	public StringComparator (Short tagName) {
		this.tagName=tagName;
	}
	
	public int compare(final Object recordOneTagValue, final Object recordTwoTagValue) {
		int cmp=0;
		String v1 = getTagValue(recordOneTagValue);
		String v2 = getTagValue(recordTwoTagValue);
		cmp=v1.compareTo(v2);
		return (cmp);
	}
	
	private String getTagValue(Object val) {
		if (val == null) return null;
        if (val instanceof String) {
            return (String)val;
        }
        throw new SAMException("Value for tag " + tagName + " is not a String: " + val.getClass());
	}
}
