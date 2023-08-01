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

import java.util.Iterator;

import org.broadinstitute.dropseqrna.utils.FilteredIterator;

import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContext;

public class AlleleFrequencyTagFilter extends FilteredIterator<VariantContext> {

	private final Log log = Log.getInstance(AlleleFrequencyTagFilter.class);
	
	private final String alleleFreqTag;
	private final double MISSING_VALUE=-1d;
	private final double minimumMinorAlleleFreq;

	public AlleleFrequencyTagFilter(final Iterator<VariantContext> underlyingIterator, final String alleleFreqTag, final Double minimumMinorAlleleFreq) {
		super(underlyingIterator);
		this.alleleFreqTag=alleleFreqTag;
		if (minimumMinorAlleleFreq==null) this.minimumMinorAlleleFreq=0; else this.minimumMinorAlleleFreq=minimumMinorAlleleFreq;
	}

	@Override
	public boolean filterOut(final VariantContext rec) {
		double minorAlleleFreq= rec.getAttributeAsDouble(this.alleleFreqTag, MISSING_VALUE);
		if (minorAlleleFreq==MISSING_VALUE) return true;
		// discard SNPs where the AF is below the limit.
		if (minorAlleleFreq < this.minimumMinorAlleleFreq) return true;
		return false;
	}

	@Override
	public void logFilterResults() {
		String msg = String.format("Minor allele freqeuncy threshold [%f] records pass [%d] records fail [%d] ",this.minimumMinorAlleleFreq, this.getRecordsPassed(), this.getRecordsFailed());  
		log.info(msg);
		
	}

}
