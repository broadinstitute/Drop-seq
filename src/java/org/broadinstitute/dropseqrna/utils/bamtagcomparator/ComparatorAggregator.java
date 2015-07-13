package org.broadinstitute.dropseqrna.utils.bamtagcomparator;

import java.util.List;

/**
 * Holds a list of one of more BAM tag comparators.
 * @author nemesh
 *
 */
public class ComparatorAggregator {

	private final List<TagComparator> comparators;
	
	private final boolean recycleComparators;
	
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
