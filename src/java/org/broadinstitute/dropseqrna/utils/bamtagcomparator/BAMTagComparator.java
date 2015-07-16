package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordComparator;
import htsjdk.samtools.SAMTagUtil;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class BAMTagComparator implements SAMRecordComparator {

	private List<Short> tagNames;

	public BAMTagComparator(Short tagName) {
		this.tagNames=new ArrayList<Short>(1);
		this.tagNames.add(tagName);
	}
	
	public BAMTagComparator (String tagName) {
		this(SAMTagUtil.getSingleton().makeBinaryTag(tagName));
	}
	
	public BAMTagComparator (Collection<Short> tags) {
		this.tagNames = new ArrayList<Short>(tags.size());
		this.tagNames.addAll(tags);
	}
	
	public BAMTagComparator (List<String> tags) {
		tagNames = new ArrayList<Short>(tagNames.size());
		for (String t: tags) {
			tagNames.add(SAMTagUtil.getSingleton().makeBinaryTag(t));
		}		
	}

	// comparator aggregater
	// list of sam record comparators.
	
	/**
	 * This only works for tags that encode Strings (for now)
	 */
	public int compare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
		int cmp=0;
		for (Short tagName: tagNames) {
			
			final String hitIndex1 = this.getTagValueAsString(tagName, samRecord1.getAttribute(tagName));
			final String hitIndex2 = this.getTagValueAsString(tagName, samRecord2.getAttribute(tagName));
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
	
	private String getTagValueAsString(Short tag, Object val) {
		if (val == null) return null;
        if (val instanceof String) {
            return (String)val;
        }
        throw new SAMException("Value for tag " + tag + " is not a String: " + val.getClass());
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
