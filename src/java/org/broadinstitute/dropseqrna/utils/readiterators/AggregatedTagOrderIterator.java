package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * A companion class to TagOrderIterator.  This class batches reads for a set of sorting tags into a batch of reads in a List.
 * @author nemesh
 *
 */
public class AggregatedTagOrderIterator implements
		CloseableIterator<Collection<SAMRecord>> {

	private final TagOrderIterator tagIter;
	
	private final List<Short> tags;

	public AggregatedTagOrderIterator(TagOrderIterator tagIter) {
		this.tagIter = tagIter;
		this.tags = tagIter.getShortTags();
	}
	
	public TagOrderIterator getTagOrderIterator () {
		return this.tagIter;
	}

	@Override
	public List<SAMRecord> next() {
		if (tagIter.hasNext() == false) return (null);
		List<SAMRecord> result = new ArrayList<SAMRecord>();

		SAMRecord r = tagIter.peek();

		List<String> currentValues = getValuesForTags(this.tags, r);
		while (tagIter.hasNext()) {
			r = tagIter.peek();

			List<String> nextValues = getValuesForTags(this.tags, r);
			if (!tagsEqual(currentValues, nextValues)) {
				break;
			}
			// this is the same set of records as before, keep going.
			// grab this record for "real" so peek gets the next record that
			// might be in the same gene.
			tagIter.next();
			// this.progress.record(r);
			result.add(r);

		}
		return (result);
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException("Remove not supported");
	}

	@Override
	public void close() {
		this.tagIter.close();
	}

	@Override
	public boolean hasNext() {
		return tagIter.hasNext();
	}

	public static List<String> getValuesForTags(List<Short> tags, SAMRecord r) {
		List<String> currentValues = new ArrayList<String>();
		for (Short t : tags) {
			Object o = r.getAttribute(t);
			if (o instanceof String)
				currentValues.add((String) o);
		}
		return (currentValues);
	}

	public static boolean tagsEqual(List<String> original, List<String> next) {
		for (int i = 0; i < original.size(); i++) {
			String s1 = original.get(i);
			String s2 = next.get(i);
			if (!s1.equals(s2))
				return (false);
		}
		return (true);
	}

}
