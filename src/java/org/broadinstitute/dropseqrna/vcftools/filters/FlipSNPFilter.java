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

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;

import java.util.Iterator;

/**
 * Filter out variants that have the potential to be flipped by a strand-flip bug.
 * These are A/T and C/G SNPs, which are rare, but vulnerable.
 * @author nemesh
 *
 */
public class FlipSNPFilter extends FilteredIterator <VariantContext> {

	public FlipSNPFilter(final Iterator<VariantContext> underlyingIterator) {
		super(underlyingIterator);
	}

	@Override
	public boolean filterOut(final VariantContext rec) {
		String ref = rec.getReference().getBaseString();
		Allele altAllele = rec.getAltAlleleWithHighestAlleleCount();
		if (altAllele==null) return false;
		String alt = altAllele.getBaseString();
		if (ref.equals("A") && alt.equals("T")) return true;
		if (ref.equals("T") && alt.equals("A")) return true;
		if (ref.equals("C") && alt.equals("G")) return true;
		if (ref.equals("G") && alt.equals("C")) return true;
		return false;
	}

}
