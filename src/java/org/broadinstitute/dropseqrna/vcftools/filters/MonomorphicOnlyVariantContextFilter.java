/*
 * MIT License
 *
 * Copyright 2018 Broad Institute
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

import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

import java.util.Iterator;
import java.util.List;

public class MonomorphicOnlyVariantContextFilter extends FilteredIterator<VariantContext> {

	private static final Log log = Log.getInstance(MonomorphicOnlyVariantContextFilter.class);
	
	final List<String> vcfSamples;

	public MonomorphicOnlyVariantContextFilter(final Iterator<VariantContext> underlyingIterator,
			final List<String> vcfSamples) {
		super(underlyingIterator);
		this.vcfSamples = vcfSamples;
	}

	@Override
	public boolean filterOut(final VariantContext rec) {
		ObjectCounter<Integer> c = new ObjectCounter<>();
		for (String s : vcfSamples) {
			Genotype g = rec.getGenotype(s);
			if (g != null) {
				if (g.isHomRef())
					c.increment(0);
				if (g.isHomVar())
					c.increment(2);
				if (g.isHet())
					c.increment(1);
			}
		}
		// filter out if not entire homRef or homVar.
		int total = c.getTotalCount();
		if (c.getCountForKey(0) == total)
			return false;
		if (c.getCountForKey(2) == total)
			return false;
		return true;
	}

	@Override
	public void logFilterResults() {
		String msg = String.format("Records pass [%d] records fail [%d] ",this.getRecordsPassed(), this.getRecordsFailed());  
		log.info(msg);						
	}

}
