package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import org.broadinstitute.dropseqrna.annotation.FunctionalData;
import org.broadinstitute.dropseqrna.annotation.FunctionalDataProcessor;
import org.broadinstitute.dropseqrna.barnyard.Utils;

import htsjdk.samtools.SAMRecord;
import picard.annotation.LocusFunction;

/**
 * This class abstracts out the logic needed by GeneFunctionIteratorWrapper, but removes the iterator
 * so that the logic can be reused on reads that are generated in different access patterns
 * @author nemesh
 *
 */
public class GeneFunctionProcessor {

	// the delimiter for BAM TAGs
	private static final String DELIMITER = ",";

	private final String geneTag;
	private final String strandTag;
	private final String functionTag;
	private final boolean assignReadsToAllGenes;
	private final FunctionalDataProcessor fdp;
	
	public GeneFunctionProcessor(final String geneTag, final String strandTag, final String functionTag,
			final boolean assignReadsToAllGenes, final StrandStrategy strandFilterStrategy,
			final Collection<LocusFunction> acceptedLociFunctions) {
		

		this.geneTag = geneTag;
		this.strandTag = strandTag;
		this.functionTag = functionTag;
		this.assignReadsToAllGenes = assignReadsToAllGenes;
		this.fdp = new FunctionalDataProcessor(strandFilterStrategy,
				acceptedLociFunctions);
	}
	
	/**
	 * Process the read to test that it passes filters, then the functional annotations are limited to a single tuple
	 * of gene/strand/function for downstream interpretation.
	 * 
	 * @param r The read to process
	 * @return The output read(s).  
	 */
	public List<SAMRecord> processRead (final SAMRecord r) {
		List<SAMRecord> result = new ArrayList<SAMRecord>();

		
		List<FunctionalData> fdList = getReadFunctions (r);
		
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
	
	public List<FunctionalData> getReadFunctions (final SAMRecord r) {
		String geneList = r.getStringAttribute(this.geneTag);
		String strandList = r.getStringAttribute(this.strandTag);
		String functionList = r.getStringAttribute(this.functionTag);

		// If you're missing the gene, you can't use this read.
		// If care about strand, and you're missing the strand, you can't use this read.
		// If care about function, and you're missing the  function, you can't use this read.
		if ((geneList == null) ||
				(fdp.getStrandStrategy() != null && strandList == null) ||
				(!fdp.getFunctions().isEmpty() && functionList == null)){
			return Collections.emptyList();
		}

		// there's at least one good copy of the read. Does the read match on
		// strand/gene, or is it assigned to multiple genes?
		final String[] genes = geneList.split(DELIMITER);
		final String[] strands = (strandList == null? null: strandList.split(DELIMITER));
		final LocusFunction[] locusFunctions = (functionList == null? null: getLocusFunctionFromRead(functionList));

		List<FunctionalData> fdList = fdp.getFilteredFunctionalData(genes,
				strands, locusFunctions, r.getReadNegativeStrandFlag());
		return (fdList);
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

	public static LocusFunction[] getLocusFunctionFromRead(final String functionList) {
		String[] fl = functionList.split(DELIMITER);
		LocusFunction[] result = new LocusFunction[fl.length];
		for (int i = 0; i < fl.length; i++) {
			LocusFunction lf = LocusFunction.valueOf(fl[i]);
			result[i] = lf;
		}
		return result;
	}
}