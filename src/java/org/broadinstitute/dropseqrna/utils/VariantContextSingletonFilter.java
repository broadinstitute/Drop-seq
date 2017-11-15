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
package org.broadinstitute.dropseqrna.utils;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.Iterator;

public class VariantContextSingletonFilter extends FilteredIterator<VariantContext>{

	private final boolean hetVarOnly;

	/**
	 * Filters out sites that have > 1 sample with an alternate allele
	 * @param underlyingIterator
	 * @param hetVarOnly True if only heterozygous sites are allowed for a variant site.  If true, observing a sample with the alternate allele will reject the record.
	 * If false then either the heterozygous or homozygous alternate will count as a variant observation.
	 * Otherwise count heterozygous and hom_var sites.
	 */
	public VariantContextSingletonFilter(final Iterator<VariantContext> underlyingIterator, final boolean hetVarOnly) {
		super(underlyingIterator);
		this.hetVarOnly=hetVarOnly;

	}

	@Override
	/**
	 * For each variant context, look at the genotypes of the samples, and count the number of samples that have
	 * an alternate allele.  If that number of samples!=1, return true to filter this record.
	 * @param rec
	 * @return
	 */
	public boolean filterOut(final VariantContext rec) {
		GenotypesContext gc = rec.getGenotypes();
		Iterator<Genotype> iter = gc.iterator();
		int count=0;
		while (iter.hasNext()) {
			Genotype g = iter.next();
			// boolean isHet = g.isHet();
			// boolean homRef = g.isHomRef();

			if (hetVarOnly & g.isHomVar())
				return true; // filter when het only and observe hom_var.
			if (g.isHet() || g.isHomVar())
				count++;
		}
		if (count==1) return false;
		return true;
	}

}
