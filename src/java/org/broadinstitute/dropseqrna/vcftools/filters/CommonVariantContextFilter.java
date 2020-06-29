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

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

import java.util.Iterator;
import java.util.List;

public class CommonVariantContextFilter extends FilteredIterator<VariantContext>{

	private final List<String> vcfSamples;
	private final int numVariantSamples;

	public CommonVariantContextFilter(final Iterator<VariantContext> underlyingIterator, final List<String> vcfSamples, final int numVariantSamples) {
		super(underlyingIterator);
		this.vcfSamples=vcfSamples;
		this.numVariantSamples=numVariantSamples;
	}

	@Override
	public boolean filterOut(final VariantContext rec) {
		ObjectCounter<Integer> c = new ObjectCounter<>();
		for (String s: vcfSamples) {
			Genotype g = rec.getGenotype(s);
			if (g!=null) {
				if (g.isHomRef()) c.increment(0);
				if (g.isHomVar()) c.increment(2);
				if (g.isHet()) c.increment(1);
			}
		}
		int sum=0;
		// if the 0 class is dominant, sum up the 1+2.
		if (c.getCountForKey(0)>=c.getCountForKey(2))
			sum = c.getCountForKey(1)+c.getCountForKey(2);
		else
			sum = c.getCountForKey(1)+c.getCountForKey(0);
		return (sum<this.numVariantSamples);
	}
}
