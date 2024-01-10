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

import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalData;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;
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
	
	private static final Log log = Log.getInstance(GeneFunctionFilteringIterator.class);
	private final GeneFunctionProcessor p;
	
	public GeneFunctionFilteringIterator(final Iterator<SAMRecord> underlyingIterator, final String geneTag, final String strandTag, 
			final String functionTag, final StrandStrategy strandFilterStrategy, final Collection<LocusFunction> acceptedLociFunctions) {
		
		super(underlyingIterator);		
		p = new GeneFunctionProcessor(geneTag, strandTag, functionTag, false, strandFilterStrategy, acceptedLociFunctions);
		
	}

	@Override
 	public boolean filterOut(final SAMRecord r) {	 							
		List<FunctionalData> fdList = p.getReadFunctions (r);
		// If there's no functional data that passes the filters, filter the read.
		if (fdList.size()==0) return true;
		// Otherwise, accept the read
		return false;
	}

	@Override
	public void logFilterResults() {
		String msg = String.format("Records pass [%d] records fail [%d] ",this.getRecordsPassed(), this.getRecordsFailed());  
		log.info(msg);										
	}
	
	
}
