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
package org.broadinstitute.dropseqrna.vcftools.filters;

import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import org.broadinstitute.dropseqrna.utils.FilteredIterator;

import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContext;

public class ChromosomeVariantFilter extends FilteredIterator<VariantContext>{

	private static final Log log = Log.getInstance(ChromosomeVariantFilter.class);
	private final Set<String> contigsToFilter;

	public ChromosomeVariantFilter(final Iterator<VariantContext> underlyingIterator, final Collection<String> contigsToFilter) {
		super(underlyingIterator);
		if (contigsToFilter!=null)
			this.contigsToFilter=new HashSet<>(contigsToFilter);
		else
			this.contigsToFilter=null;
	}

	@Override
	public boolean filterOut(final VariantContext rec) {
		// short circuit if there are no contigs to filter.
		if (this.contigsToFilter==null) return false;
		return this.contigsToFilter.contains(rec.getContig());
	}

	@Override
	public void logFilterResults() {
		String msg = String.format("Records pass [%d] records fail [%d] ",this.getRecordsPassed(), this.getRecordsFailed());  
		log.info(msg);										
	}

}
