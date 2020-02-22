package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;
import java.util.function.Predicate;

import htsjdk.samtools.SAMRecord;

/**
 * Given a tag name and a collection of values, filters reads that do not have one of the values. 
 * @author nemesh
 *
 */
public class RequiredTagStringValuePredicate implements Predicate<SAMRecord> {
	private final String requiredTag;
	private final Set<String> values;
	private final boolean exclude;
	
    public RequiredTagStringValuePredicate(final String requiredTag, final Collection<String> tagValues, boolean exclude) {
    	if (tagValues==null) values=null;
    	else this.values=new HashSet<>(tagValues);
    	this.requiredTag=requiredTag;
        this.exclude=exclude;
    }

    @Override
    public boolean test(SAMRecord rec) {
    	// if values are null or empty, don't filter.
    	if (values==null) return false;
		if (values.isEmpty()) return false;
		
		boolean contains=values.contains(rec.getStringAttribute(this.requiredTag));
		// if exclude filter if contained.
    	if (exclude) return !contains;
    	// if include, filter if not contained 
    	return contains;    	    	
    }
		
}
