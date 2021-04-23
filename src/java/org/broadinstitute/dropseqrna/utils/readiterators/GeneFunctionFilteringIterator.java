/*
 * MIT License
 *
 * Copyright 2021 Broad Institute
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
import java.util.Iterator;
import java.util.List;

import org.broadinstitute.dropseqrna.annotation.FunctionalData;
import org.broadinstitute.dropseqrna.annotation.FunctionalDataProcessor;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;

import htsjdk.samtools.SAMRecord;
import picard.annotation.LocusFunction;

/**
 * This is a simplified version of the logic used for GeneFunctionIteratorWrapper
 * However, this does not retag the read with a simplified preferred tag.
 * Instead, if the read has the proper functional information to be considered for downstream analysis, keep the read as-is
 * Otherwise, filter the read.
 * @author nemesh
 *
 */
public class GeneFunctionFilteringIterator extends FilteredIterator<SAMRecord> {
	
	private static final String DELIMITER = ",";
	private final String geneTag;
	private final String strandTag;
	private final String functionTag;	
	private final FunctionalDataProcessor fdp;
	
	public GeneFunctionFilteringIterator(final Iterator<SAMRecord> underlyingIterator, final String geneTag, final String strandTag, 
			final String functionTag, final StrandStrategy strandFilterStrategy, final Collection<LocusFunction> acceptedLociFunctions) {
		
		super(underlyingIterator);
		this.geneTag = geneTag;
		this.strandTag = strandTag;
		this.functionTag = functionTag;
		this.fdp = new FunctionalDataProcessor(strandFilterStrategy,
				acceptedLociFunctions);
	}

	
	@Override
 	public boolean filterOut(final SAMRecord r) {
	 	String geneList = r.getStringAttribute(this.geneTag);
		String strandList = r.getStringAttribute(this.strandTag);
		String functionList = r.getStringAttribute(this.functionTag);

		// If you're missing the gene, you can't use this read.
		// If care about strand, and you're missing the strand, you can't use this read.
		// If care about function, and you're missing the  function, you can't use this read.
		if ((geneList == null) ||
				(fdp.getStrandStrategy() != null && strandList == null) ||
				(!fdp.getFunctions().isEmpty() && functionList == null)){
			return true;
		}
		
		// there's at least one good copy of the read. Does the read match on
		// strand/gene, or is it assigned to multiple genes?
		final String[] genes = geneList.split(DELIMITER);
		final String[] strands = (strandList == null? null: strandList.split(DELIMITER));
		final LocusFunction[] locusFunctions = (functionList == null? null: GeneFunctionIteratorWrapper.getLocusFunctionFromRead(functionList));

		List<FunctionalData> fdList = fdp.getFilteredFunctionalData(genes,
				strands, locusFunctions, r.getReadNegativeStrandFlag());

		// If there's no functional data that passes the filters, filter the read.
		if (fdList.size()==0) return true;
		// Otherwise, accept the read
		return false;
	}
}
