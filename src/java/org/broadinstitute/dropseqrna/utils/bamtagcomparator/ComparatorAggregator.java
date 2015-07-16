package org.broadinstitute.dropseqrna.utils.bamtagcomparator;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Holds a list of one of more BAM tag comparators.
 * These define how the tags are sorted.
 * There are two options for this list of comparators
 * 1) List one or more comparators, and recycle comparators - if you have 2 defined comparators and 3 tags, then tag 3 gets the first comparator as they are recycled
 * 2) recycleComparators is false, and you have exactly 1 comparator per SNP tag.
 * @author nemesh
 *
 */
public class ComparatorAggregator {

	private final List<TagComparator> comparators;
	
	private final boolean recycleComparators;
	
	public ComparatorAggregator(TagComparator comparator, boolean recycleComparators) {
		List<TagComparator> temp = new ArrayList<TagComparator>();
		temp.add(comparator);
		this.comparators=temp;
		this.recycleComparators=recycleComparators;
	}
	
	public ComparatorAggregator (TagComparator [] comparators, boolean recycleComparators) {
		this.comparators=Arrays.asList(comparators);
		this.recycleComparators=recycleComparators;
	}
	
	public ComparatorAggregator (List<TagComparator> comparators, boolean recycleComparators) {
		this.recycleComparators=recycleComparators;
		this.comparators=comparators;
	}
	
	public List<TagComparator> getComparators() {
		return this.comparators;
	}
	
	public TagComparator get(int index) {
		if (recycleComparators) {
			index=index%comparators.size();
		}
		return this.comparators.get(index);
	}
	
}
