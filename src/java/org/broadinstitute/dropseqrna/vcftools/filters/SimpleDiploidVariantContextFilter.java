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

import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.broadinstitute.dropseqrna.utils.FilteredIterator;

import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContext.Type;

public class SimpleDiploidVariantContextFilter extends FilteredIterator <VariantContext> {

	private final Log log = Log.getInstance(SimpleDiploidVariantContextFilter.class);

	private final boolean filterNonSNPs;
	private final boolean filterFilterFlagedVariants;
	private final int maxNumAlleles;
	private final boolean retainMonmorphicSNPs;
	private final Integer maxAlleleLength;
	private final boolean verbose;
	private final Set<String> canonicalBaseSet= new HashSet<>(Arrays.asList("A", "C", "G", "T"));
	private final List<byte []> canonicalBaseList;
	
	/**
	 * Filters the base length of alleles to 1 by default.
	 * @param underlyingIterator
	 * @param filterNonSNPs
	 * @param filterFilterFlagedVariants
	 * @param maxNumAlleles
	 * @param retainMonmorphicSNPs
	 */
	public SimpleDiploidVariantContextFilter (final Iterator<VariantContext> underlyingIterator, final boolean filterNonSNPs, final boolean filterFilterFlagedVariants, final int maxNumAlleles, final boolean retainMonmorphicSNPs) {
		this(underlyingIterator, filterNonSNPs, filterFilterFlagedVariants, maxNumAlleles, retainMonmorphicSNPs, 1, false);
	}

	/**
	 *
	 * @param underlyingIterator
	 * @param filterNonSNPs
	 * @param filterFilterFlagedVariants
	 * @param maxNumAlleles
	 * @param retainMonmorphicSNPs
	 * @param maxAlleleLength if set to null, do not filter on allele length.
	 */
	public SimpleDiploidVariantContextFilter (final Iterator<VariantContext> underlyingIterator, final boolean filterNonSNPs, final boolean filterFilterFlagedVariants, final int maxNumAlleles, final boolean retainMonmorphicSNPs, final Integer maxAlleleLength, final boolean verbose) {
		super(underlyingIterator);
		this.filterNonSNPs=filterNonSNPs;
		this.filterFilterFlagedVariants=filterFilterFlagedVariants;
		this.maxNumAlleles=maxNumAlleles;
		this.retainMonmorphicSNPs=retainMonmorphicSNPs;
		this.maxAlleleLength=maxAlleleLength;
		this.verbose=verbose;
		// construct the canonical bases as bytes [].
		canonicalBaseList = canonicalBaseSet.stream().map(x->x.toUpperCase().getBytes()).collect(Collectors.toList());	
		
	}
	public SimpleDiploidVariantContextFilter (final Iterator<VariantContext> underlyingIterator) {
		this(underlyingIterator, true, true, 2, false);
	}

	@Override
	public boolean filterOut(final VariantContext site) {
		// if requested, filter out any "filtered" site.
		if (filterFilterFlagedVariants && site.isFiltered()) {
			if (verbose) log.info("Rejecting variant site filtered "+site.toStringWithoutGenotypes());
			return true;

		}
		// if requested, filter SNPs with greater than <maxNumAlleles> alleles.
		List<Allele> alleles= site.getAlleles();
		for (Allele a: alleles)
			if (this.maxAlleleLength!=null && a.length()>this.maxAlleleLength) {
				if (verbose) log.info("Rejecting variant alleles too long "+site.toStringWithoutGenotypes());
				return true;
			}
		if (site.getAlleles().size()>this.maxNumAlleles) {
			if (verbose) log.info("Rejecting variant too many alleles "+site.toStringWithoutGenotypes());
			return true;
		}

		// look for non-canonical bases: One of the alleles is "*" or "N".
		for (Allele a: alleles) {
			if (!isCanonicalAllele(a)) {
				if (verbose) log.info("Rejecting variant because allele is non-canonical "+site.toStringWithoutGenotypes());
				return true;
			}
		}
		
		// if requested, filter out nonSNP sites.
		// if the site is a non snp because it's monomorphic and you want to retain those, then return false.
		if (filterNonSNPs && !site.isSNP()) {
			if (site.getType()==Type.NO_VARIATION && retainMonmorphicSNPs)
				return false;
			if (verbose) log.info("Rejecting variant not a SNP or monomorphic in population "+site.toStringWithoutGenotypes());
			return true;
		}

		return false;
	}
	
	private boolean isCanonicalAllele (Allele a) {
		for (byte [] b: canonicalBaseList) {
			if (a.basesMatch(b))
				return true;
		}
		return false;
		
	}
}


