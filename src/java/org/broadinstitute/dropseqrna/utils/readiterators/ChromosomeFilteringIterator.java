package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.broadinstitute.dropseqrna.utils.FilteredIterator;

import htsjdk.samtools.SAMRecord;

public class ChromosomeFilteringIterator extends FilteredIterator<SAMRecord>{

	private final Set<String> contigsToFilter;
	private final boolean excludeContig;

	/**
	 * Filter out records that have a chromosome that is contained in the contigsToFilter collection.
	 * @param underlyingIterator
	 * @param contigsToFilter
	 */
	public ChromosomeFilteringIterator(final Iterator<SAMRecord> underlyingIterator, final Collection<String> contigsToFilter) {
		this(underlyingIterator, contigsToFilter, true);
	}

	/**
	 * Exclude or include chromosomes.
	 *
	 * @param underlyingIterator The iterator to filter
	 * @param contigsToFilter The list of contigs to filter
	 * @param excludeContig If true, then filter out all records that contain a contig in the contigsToFilter set.  If false, filter out all records not contain in the contigsToFilter set.
	 */
	public ChromosomeFilteringIterator(final Iterator<SAMRecord> underlyingIterator, final Collection<String> contigsToFilter, final boolean excludeContig) {
		super(underlyingIterator);
		if (contigsToFilter!=null)
			this.contigsToFilter=new HashSet<>(contigsToFilter);
		else
			this.contigsToFilter=null;
		this.excludeContig=excludeContig;
	}


	@Override
	public boolean filterOut(final SAMRecord rec) {
		// short circuit if there are no contigs to filter.
		if (this.contigsToFilter==null) return false;
		// if you're excluding contigs, then return true to filter if the record is contained in the contigs
		if (this.excludeContig) return this.contigsToFilter.contains(rec.getContig());
		// if you're including contigs, then return false if the record is contained in the contigs.
		return !this.contigsToFilter.contains(rec.getContig());

	}


}
