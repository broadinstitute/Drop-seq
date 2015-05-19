package org.broadinstitute.dropseqrna.barnyard;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.ProgressLogger;

import java.util.ArrayList;
import java.util.List;

public class Utils {
	private static Log log = Log.getInstance(Utils.class);
	private ProgressLogger progress = new ProgressLogger(log, 1000000);
	
	
	
	
	/**
	 * When you sort an iterator on multiple fields and some of the fields [1 or more] may be null, this skips past all the null entries and starts with the first non-null entry.
	 * @param iter The iterator to...iterate on
	 * @param emptyAttribute The attribute(s) of the read to scan.
	 * 
	 * @return
	 */
	public PeekableIterator<SAMRecord> primeIterator (PeekableIterator<SAMRecord> iter, String...emptyAttribute) {
		ProgressLogger primeLog = new ProgressLogger(log, 1000000, "Skipped records without tags "+ getFormattedString(emptyAttribute));
		if (iter.hasNext()==false) return (iter);
		
		SAMRecord r = iter.peek();
		// seek to the first gene.
		while (iter.hasNext()) {
			r=iter.peek();
			int numNotNull=0;
			for (String key: emptyAttribute) {
				Object value = r.getAttribute(key);
				if (value!=null) numNotNull++;
			}
			if (numNotNull==emptyAttribute.length){ 
				break;
			}
			r=iter.next();
			primeLog.record(r);
		}
		return (iter);
	}
	
	private String getFormattedString (String...x) {
		StringBuilder b= new StringBuilder();
		for (int i=0; i<x.length; i++) {
			b.append(x[i]);
			if (i<x.length) {
				b.append(",");
			}
		}
		return (b.toString());
	}
	
	
	/**
	 * Get all the reads for one set of tags on an iterator.
	 * When this method finishes, it should have peeked the next read that has a different tag than the tags listed.
	 * @param iter A peekable iterator of reads in gene tag order.
	 * @return A list of records, or null if the iterator is out.  Returns an empty list if no reads satisfy the filter critieria for this batch.
	 */
	public List<SAMRecord> getReadsForNextBatch (PeekableIterator<SAMRecord> iter, List<String> tags, List<String> filterTagsSet, Integer mapQuality) {
		if (iter.hasNext()==false) return (null);
		List<SAMRecord> result = new ArrayList<SAMRecord>();
		
		SAMRecord r = iter.peek();
		
		
		List<String> currentValues = getValuesForTags(tags, r);
		while (iter.hasNext()) {
			r=iter.peek();
			
			List<String> nextValues = getValuesForTags(tags, r);
			if (testTagsNotEqual(currentValues, nextValues)) {
				break;
			}
			// this is the same set of records as before, keep going.
			// grab this record for "real" so peek gets the next record that might be in the same gene.
			iter.next();
			this.progress.record(r);
			
			boolean tagsSet=testTagSet(filterTagsSet, r);
			// filter read to be a primary read of the right map quality that is exonic, or don't add it.
			if (!r.isSecondaryOrSupplementary()&& r.getMappingQuality()>=mapQuality && tagsSet) {
				result.add(r);
			}
		}
		return (result);
		
		
	}
	
	private boolean testTagSet (List<String> tags, SAMRecord r) {
		if (tags==null || tags.isEmpty()) return true;
		for (String t: tags) {
			Object o = r.getAttribute(t);
			if (o==null) return (false);
		}
		return (true);
		
	}
	
	private List<String> getValuesForTags(List<String>tags, SAMRecord r) {
		List<String> currentValues = new ArrayList<String>();
		for (String t: tags) {
			currentValues.add(r.getStringAttribute(t));
		}
		return (currentValues);
	}
	
	private boolean testTagsNotEqual (List<String> original, List<String> next) {
		for (int i=0; i<original.size(); i++) {
			String s1 = original.get(i);
			String s2 = next.get(i);
			if (!s1.equals(s2)) return (true);
		}
		return (false);
	}
	
	
}
