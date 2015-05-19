package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordComparator;

import java.util.ArrayList;
import java.util.List;

public class BAMTagComparator implements SAMRecordComparator {

	private List<String> tagNames;

	public BAMTagComparator(String tagName) {
		this.tagNames=new ArrayList<String>(1);
		this.tagNames.add(tagName);
	}
	
	public BAMTagComparator(List<String> tagNames) {
		this.tagNames=tagNames;
	}

	public int compare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
		int cmp=0;
		for (String tagName: tagNames) {
			final String hitIndex1 = samRecord1.getStringAttribute(tagName);
			final String hitIndex2 = samRecord2.getStringAttribute(tagName);
			if (hitIndex1 != null) {
				if (hitIndex2 == null)
					return 1;
				else {
					cmp = hitIndex1.compareTo(hitIndex2);
					if (cmp != 0)
						return cmp;
				}
			} else if (hitIndex2 != null) {
				return -1;
			}	
		}
		return 0;
	}

	/**
	 * Less stringent compare method than the regular compare. If the two
	 * records are equal enough that their ordering in a sorted SAM file would
	 * be arbitrary, this method returns 0.
	 * 
	 * @return negative if samRecord1 < samRecord2, 0 if equal, else positive
	 */
	public int fileOrderCompare(final SAMRecord samRecord1,
			final SAMRecord samRecord2) {
		return (0);
	}

	
}
