package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTagUtil;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

/**
 * Keeps or rejects a read if the set tag has a value in the set of potential values.
 * @author nemesh
 *
 */
public class TagValueProcessor implements SAMReadProcessorI {

	private final Short tag;
	private final Set<String> values;
	private boolean keepReadWithValue;
	
	/**
	 * Set up the processor with a BAM tag to filter reads on.
	 * 
	 * @param tag The BAM tag to check String values of.
	 * @param values The collection of String values that will be accepted or rejected.
	 * @param keepReadWithValue
	 */
	public TagValueProcessor (String tag, Collection<String> values, boolean keepReadWithValue) {
		this.tag = SAMTagUtil.getSingleton().makeBinaryTag(tag);
		this.values = new HashSet<String>(values);
		this.keepReadWithValue = keepReadWithValue;
	}
	
	@Override
	public Collection<SAMRecord> processRead(SAMRecord r,
			Collection<SAMRecord> outList) {
		
		Object o = r.getAttribute(this.tag);
		if (o instanceof String) {
			String value = (String) o;
			if (values.contains(value)==keepReadWithValue) {
				outList.add(r);
			}
		}		
		return outList;
	}

}
