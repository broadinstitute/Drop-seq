package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.Collection;
import java.util.HashSet;
import java.util.Set;
import java.util.function.Predicate;

import htsjdk.samtools.SAMRecord;

public class ChromosomeFilteringPredicate extends RequiredTagPredicate implements Predicate<SAMRecord> {

	
	private final Set<String> contigsToFilter;
	private final boolean excludeContig;
	
	public ChromosomeFilteringPredicate(final Collection<String> contigsToFilter, final boolean excludeContig) {
		if (contigsToFilter!=null)
			this.contigsToFilter=new HashSet<>(contigsToFilter);
		else
			this.contigsToFilter=null;
		this.excludeContig=excludeContig;		
	}
	
	@Override
    public boolean test(SAMRecord rec) {
		// short circuit if there are no contigs to filter.
		if (this.contigsToFilter==null) return true;
		if (this.contigsToFilter.isEmpty()) return true;
		
		// if you're excluding contigs, then return true to filter if the record is contained in the contigs
		String contig = rec.getContig();
		if (this.excludeContig) return !this.contigsToFilter.contains(contig);
		// if you're including contigs, then return true if the record is contained in the contigs.
		return this.contigsToFilter.contains(contig);
	}
	
}
