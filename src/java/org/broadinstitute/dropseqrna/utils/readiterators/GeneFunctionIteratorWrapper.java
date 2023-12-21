package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.dropseqrna.utils.CountChangingIteratorWrapper;

import picard.annotation.LocusFunction;

import java.util.*;

public class GeneFunctionIteratorWrapper extends
		CountChangingIteratorWrapper<SAMRecord>  {



	private final GeneFunctionProcessor geneFunctionProcessor;
	
	public GeneFunctionIteratorWrapper(
			final Iterator<SAMRecord> underlyingIterator, final String geneTag,
			final String strandTag, final String functionTag,
			final boolean assignReadsToAllGenes,
			final StrandStrategy strandFilterStrategy,
			final Collection<LocusFunction> acceptedLociFunctions) {
		super(underlyingIterator);
		geneFunctionProcessor = new GeneFunctionProcessor(geneTag, strandTag, functionTag, assignReadsToAllGenes, strandFilterStrategy, acceptedLociFunctions);
	}

	@Override
	public void processRecord(final SAMRecord r) {
		List<SAMRecord> result = geneFunctionProcessor.processRead(r);
		if (result.size()==0)			
			return;
		
		for (SAMRecord rr: result)
			queueRecordForOutput(rr);
	}

	public GeneFunctionProcessor getGeneFunctionProcessor() {
		return geneFunctionProcessor;
	}


}
