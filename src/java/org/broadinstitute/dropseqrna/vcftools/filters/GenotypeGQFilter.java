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
import htsjdk.variant.variantcontext.Genotype;

/**
 * Filters Genotype objects by a given genotype quality.
 * Optionally filter out no call genotypes.
 * @author nemesh
 *
 */
public class GenotypeGQFilter extends FilteredIterator <Genotype> {

	// private static final Log log = Log.getInstance(GenotypeGQFilter.class);
	
	private final int GQThreshold;
	private boolean filterUncalled;

	public GenotypeGQFilter (final Iterator<Genotype> underlyingIterator, final int genotypeQualityThreshold) {
		this (underlyingIterator, genotypeQualityThreshold, true);
	}

	public GenotypeGQFilter (final Iterator<Genotype> underlyingIterator, final int genotypeQualityThreshold, final boolean filterUncalled) {
		super(underlyingIterator);
		this.GQThreshold=genotypeQualityThreshold;
		this.filterUncalled=filterUncalled;
	}

	@Override
	public boolean filterOut(final Genotype rec) {
		if (!rec.isCalled() & filterUncalled) return true;
		boolean flag = (rec.getGQ()<GQThreshold);
		return flag;
	}

	@Override
	public void logFilterResults() {
		// This is disabled as it is far too verbose - it's typically used per site.
		// String msg = String.format("GQ Threshold [%d] filter uncalled [%B] records pass [%d] records fail [%d] ",this.GQThreshold, this.filterUncalled, this.getRecordsPassed(), this.getRecordsFailed());  
		// log.info(msg);								
		
	}
}
