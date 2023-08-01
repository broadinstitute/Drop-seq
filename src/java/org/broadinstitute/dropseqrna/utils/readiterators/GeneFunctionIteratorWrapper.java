package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.dropseqrna.annotation.FunctionalData;
import org.broadinstitute.dropseqrna.annotation.FunctionalDataProcessor;
import org.broadinstitute.dropseqrna.barnyard.Utils;
import org.broadinstitute.dropseqrna.utils.CountChangingIteratorWrapper;
import org.broadinstitute.dropseqrna.utils.PassFailTrackingIteratorI;

import picard.annotation.LocusFunction;

import java.util.*;

public class GeneFunctionIteratorWrapper extends
		CountChangingIteratorWrapper<SAMRecord>  {

	
	private final GeneFunctionProcessor p;
	
	public GeneFunctionIteratorWrapper(
			final Iterator<SAMRecord> underlyingIterator, final String geneTag,
			final String strandTag, final String functionTag,
			final boolean assignReadsToAllGenes,
			final StrandStrategy strandFilterStrategy,
			final Collection<LocusFunction> acceptedLociFunctions) {
		super(underlyingIterator);
		p = new GeneFunctionProcessor(geneTag, strandTag, functionTag, assignReadsToAllGenes, strandFilterStrategy, acceptedLociFunctions);
	}

	@Override
	public void processRecord(final SAMRecord r) {
		List<SAMRecord> result = p.processRead(r);
		if (result.size()==0)			
			return;
		
		for (SAMRecord rr: result)
			queueRecordForOutput(rr);
	}


	
}
