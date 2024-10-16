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

import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContext.Type;
import org.broadinstitute.dropseqrna.TranscriptomeException;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;

import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

public class SimpleDiploidVariantContextFilter extends FilteredIterator <VariantContext> {

	private final Log log = Log.getInstance(SimpleDiploidVariantContextFilter.class);

	private final boolean filterNonSNPs;
	private final boolean filterFilterFlaggedVariants;
	private final int maxNumAlleles;
	// Internally stored as int, converts null to -1.
	private final int maxAlleleLength;
	private final boolean verbose;
	private final List<Character> canonicalBaseList= Arrays.asList('A', 'C', 'G', 'T');
	private final byte [] canonicalBaseArray= new byte [canonicalBaseList.size()];

	/**
	 * Filters the base length of alleles to 1 by default.
	 * @param underlyingIterator
	 * @param filterNonSNPs
	 * @param filterFilterFlaggedVariants
	 * @param maxNumAlleles
	 */
	public SimpleDiploidVariantContextFilter (final Iterator<VariantContext> underlyingIterator, final boolean filterNonSNPs, final boolean filterFilterFlaggedVariants,
			final int maxNumAlleles) {
		this(underlyingIterator, filterNonSNPs, filterFilterFlaggedVariants, maxNumAlleles, 1, false);
	}

	/**
	 * Filter variant context sites to the subset that meet the requirements.
	 * 
	 * @param underlyingIterator An iterator of variant context objects to filter
	 * @param filterNonSNPs Remove non-SNP variants.
	 * @param filterFlaggedVariants Remove variants that dont have a PASS filter
	 * @param maxNumAlleles Restrict the maximum number of alleles.  Typically used to remove multi-allelic SNPs
	 * @param maxAlleleLength if set to null, do not filter on allele length.  If set to length one, also checks that the base for 
	 * each allele is one of the canonical bases - A/C/G/T.  This excludes N and * bases.
	 */
	public SimpleDiploidVariantContextFilter (final Iterator<VariantContext> underlyingIterator, final boolean filterNonSNPs, final boolean filterFlaggedVariants,
			final Integer maxNumAlleles, final Integer maxAlleleLength, final boolean verbose) {
		super(underlyingIterator);
		this.filterNonSNPs=filterNonSNPs;
		this.filterFilterFlaggedVariants =filterFlaggedVariants;
		
		if (maxAlleleLength==null) {
			this.maxAlleleLength=-1;
		} else {
			this.maxAlleleLength=maxAlleleLength;
		}
		
		this.maxNumAlleles=maxNumAlleles;
		this.verbose=verbose;
		
		for (int i=0; i<canonicalBaseList.size(); i++) {			
			canonicalBaseArray[i]=(byte) canonicalBaseList.get(i).charValue();
		}
		// sorted for later binary search
		Arrays.sort(canonicalBaseArray);
		
	}
	public SimpleDiploidVariantContextFilter (final Iterator<VariantContext> underlyingIterator) {
		this(underlyingIterator, true, true, 2);
	}

	@Override
	public boolean filterOut(final VariantContext site) {
		// Remove variants with symbolic alleles.
		// These would not pass the isCanonicalAllele test later.
		if (site.hasSymbolicAlleles()) {
			return true;
		}
		// if requested, filter out any "filtered" site.
		if (filterFilterFlaggedVariants && site.isFiltered()) {
			if (verbose) log.info("Rejecting variant site filtered "+site.toStringWithoutGenotypes());
			return true;
		}
		// if requested, filter SNPs with greater than <maxNumAlleles> alleles.
		List<Allele> alleles= site.getAlleles();
		// validate proper number of alleles
		if (site.getAlleles().size()>this.maxNumAlleles) {
			if (verbose) log.info("Rejecting variant too many alleles "+site.toStringWithoutGenotypes());
			return true;
		}
		
		// look for non-canonical bases: One of the alleles is "*" or "N".  Reject alleles that are longer than length one.
		for (Allele a: alleles) {
			if (!isCanonicalAllele(a)) {
				if (verbose) log.info("Rejecting variant [" +site.toStringWithoutGenotypes()+ "] because allele is non-canonical ["+a.toString()+"]");
				return true;
			}
		}
		
		// if requested, filter out nonSNP sites.
		if (filterNonSNPs && !site.isSNP()) {
			if (verbose) log.info("Rejected variant not a SNP and nonSNPs are filtered" +site.toStringWithoutGenotypes());
			return true;
		}

		// if the site does not vary and we don't want to retain monomorphic sites then filter this site.

		// VariantContext passes all filters
		return false;
	}
	
	/**
	 * Test if the allele is canonical
	 * If the max size length is not set, then don't perform any tests
	 * If the max size is set, then enforce that restriction
	 * If the max size is set to 1, then also enforce that the allele must be A,C,G,T.
	 * @param a The allele to test
	 * @return True if the allele passes these tests.
	 */
	private boolean isCanonicalAllele (Allele a) {
		// if the max allele length is null, then there's no restriction on allele length or content.
		if (this.maxAlleleLength==-1) 
			return true;

		// filter on allele length
		byte [] bases = a.getBases();
		// reject alleles with length > maxAlleleLength
		if (bases.length>maxAlleleLength)
			return false;
		
		// If allele must be of length 1 check that bases are expected canonical.
		if (maxAlleleLength==1) {		
			byte thisBase = a.getBases()[0];
			int x = Arrays.binarySearch(canonicalBaseArray, thisBase);
			return (x>=0);	
		}
		throw new TranscriptomeException("Canonical Allele filter undefined behavior, please add proper tests here.");
		
				
	}

	@Override
	public void logFilterResults() {
		String msg = String.format("Records pass [%d] records fail [%d] ",this.getRecordsPassed(), this.getRecordsFailed());  
		log.info(msg);		
	}
	
}


