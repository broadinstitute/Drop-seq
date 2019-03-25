package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.dropseqrna.annotation.FunctionalData;
import org.broadinstitute.dropseqrna.annotation.FunctionalDataProcessor;
import org.broadinstitute.dropseqrna.barnyard.Utils;
import org.broadinstitute.dropseqrna.utils.CountChangingIteratorWrapper;
import picard.annotation.LocusFunction;

import java.util.*;

public class GeneFunctionIteratorWrapper extends
		CountChangingIteratorWrapper<SAMRecord> {

	// the delimiter for BAM TAGs
	private static final String DELIMITER = ",";
	private final String geneTag;
	private final String strandTag;
	private final String functionTag;
	private final boolean assignReadsToAllGenes;
	private final FunctionalDataProcessor fdp;


	public GeneFunctionIteratorWrapper(
			final Iterator<SAMRecord> underlyingIterator, final String geneTag,
			final String strandTag, final String functionTag,
			final boolean assignReadsToAllGenes,
			final StrandStrategy strandFilterStrategy,
			final Collection<LocusFunction> acceptedLociFunctions) {
		super(underlyingIterator);

		this.geneTag = geneTag;
		this.strandTag = strandTag;
		this.functionTag = functionTag;
		this.assignReadsToAllGenes = assignReadsToAllGenes;
		this.fdp = new FunctionalDataProcessor(strandFilterStrategy,
				acceptedLociFunctions);
	}

	@Override
	public void processRecord(final SAMRecord r) {
		List<SAMRecord> result = processRead(r);
		for (SAMRecord rr: result)
			queueRecordForOutput(rr);

	}

	public List<SAMRecord> processRead (final SAMRecord r) {
		List<SAMRecord> result = new ArrayList<SAMRecord>();

		String geneList = r.getStringAttribute(this.geneTag);
		String strandList = r.getStringAttribute(this.strandTag);
		String functionList = r.getStringAttribute(this.functionTag);

		// If you're missing the gene, you can't use this read.
		// If care about strand, and you're missing the strand, you can't use this read.
		// If care about function, and you're missing the  function, you can't use this read.
		if ((geneList == null) ||
				(fdp.getStrandStrategy() != null && strandList == null) ||
				(!fdp.getFunctions().isEmpty() && functionList == null)){
			return result;
		}

		// there's at least one good copy of the read. Does the read match on
		// strand/gene, or is it assigned to multiple genes?
		final String[] genes = geneList.split(DELIMITER);
		final String[] strands = (strandList == null? null: strandList.split(DELIMITER));
		final LocusFunction[] locusFunctions = (functionList == null? null: getLocusFunctionFromRead(functionList));

		List<FunctionalData> fdList = fdp.getFilteredFunctionalData(genes,
				strands, locusFunctions, r.getReadNegativeStrandFlag());

		// If there's no functional data that passes the filters, return.
		if (fdList.size()==0) return result;

		// filter to only the preferred class to resolve reads that overlap a coding region on one gene
		// and an intronic region on a different gene.
		// this filters the data to just the coding gene.
		fdList = fdp.filterToPreferredAnnotations(fdList);

		// If there's only one functional data result, re-tag the read, and
		// queue it, then short circuit out.
		if (fdList.size() == 1) {
			FunctionalData fd = fdList.get(0);
			SAMRecord rr = assignTagsToRead(r, fd);
			result.add(rr);
			return result;
		}
		// more than 1 read.
		if (this.assignReadsToAllGenes)
			// if fdList is empty, no records are added
			for (FunctionalData fd : fdList) {
				SAMRecord rr = assignTagsToRead(Utils.getClone(r), fd);
				result.add(rr);
			}
		return result;
	}

	private SAMRecord assignTagsToRead(final SAMRecord r,
			final FunctionalData fd) {
		r.setAttribute(geneTag, fd.getGene());
		r.setAttribute(strandTag, fd.getStrand());
		if (fd.getLocusFunction() != null) {
			r.setAttribute(functionTag, fd.getLocusFunction().name());
		}
		return (r);
	}

	private LocusFunction[] getLocusFunctionFromRead(final String functionList) {
		String[] fl = functionList.split(DELIMITER);
		LocusFunction[] result = new LocusFunction[fl.length];
		for (int i = 0; i < fl.length; i++) {
			LocusFunction lf = LocusFunction.valueOf(fl[i]);
			result[i] = lf;
		}
		return result;
	}

}
