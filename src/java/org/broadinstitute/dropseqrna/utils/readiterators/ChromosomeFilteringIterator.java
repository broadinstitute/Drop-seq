/*
 * MIT License
 *
 * Copyright 2017 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
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
