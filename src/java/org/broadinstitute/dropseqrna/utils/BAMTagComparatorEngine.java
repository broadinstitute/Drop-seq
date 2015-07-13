package org.broadinstitute.dropseqrna.utils;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordComparator;
import htsjdk.samtools.SAMTagUtil;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.broadinstitute.dropseqrna.utils.bamtagcomparator.ComparatorAggregator;
import org.broadinstitute.dropseqrna.utils.bamtagcomparator.TagComparator;


public class BAMTagComparatorEngine {

	public class BAMTagComparator implements SAMRecordComparator {

		private final List<Short> tagNames;
		private final ComparatorAggregator comparatorAggregator;
		
		public BAMTagComparator(Short tagName, ComparatorAggregator comparatorAggregator) {
			this.tagNames=new ArrayList<Short>(1);
			this.tagNames.add(tagName);
			this.comparatorAggregator=comparatorAggregator;
		}
		
		public BAMTagComparator (String tagName, ComparatorAggregator comparatorAggregator) {
			this(SAMTagUtil.getSingleton().makeBinaryTag(tagName), comparatorAggregator);
		}
		
		public BAMTagComparator (Collection<Short> tags, ComparatorAggregator comparatorAggregator) {
			this.tagNames = new ArrayList<Short>(tags.size());
			this.tagNames.addAll(tags);
			this.comparatorAggregator=comparatorAggregator;
		}
		
		public BAMTagComparator (List<String> tags, ComparatorAggregator comparatorAggregator) {
			tagNames = new ArrayList<Short>(tags.size());
			for (String t: tags) {
				tagNames.add(SAMTagUtil.getSingleton().makeBinaryTag(t));
			}		
			this.comparatorAggregator=comparatorAggregator;
		}

		// comparator aggregater
		// list of sam record comparators.
		
		/**
		 * This only works for tags that encode Strings (for now)
		 */
		public int compare(final SAMRecord samRecord1, final SAMRecord samRecord2) {
			int cmp=0;
			for (int i=0; i<tagNames.size(); i++) {
				Short tagName=tagNames.get(i);
				TagComparator tc = this.comparatorAggregator.get(i);
			
				final Object o1 = samRecord1.getAttribute(tagName);
				final Object o2 = samRecord2.getAttribute(tagName);
								
				if (o1 != null) {
					if (o2 == null)
						return 1;
					else {
						cmp = tc.compare(o1, o2);
						if (cmp != 0)
							return cmp;
					}
				} else if (o2 != null) {
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

}
