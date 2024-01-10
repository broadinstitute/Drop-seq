package org.broadinstitute.dropseqrna.annotation.functionaldata;

import org.broadinstitute.dropseqrna.annotation.AnnotationUtils;
import org.broadinstitute.dropseqrna.barnyard.Utils;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import picard.annotation.LocusFunction;

import java.util.*;

/**
 * Given lists of gene names, strands of the gene, and functional annotations of the genes,
 * decide which genes/strands/functions should be retained.
 *
 * The abstract class supports minimal functions for multiple strategy / implementations
 * @author nemesh
 *
 */
public abstract class AbstractFunctionalDataProcessor {

	private final StrandStrategy strandStrategy;
	private final Collection <LocusFunction> acceptedFunctions;

	/**
	 * Initialize a processor with a set of locus functions that are acceptable, and a strand strategy.
	 *
	 * @param strandStrategy    The strategy for filtering on read and gene strand:
	 *                          StrandStrategy.SENSE retains all annotations where the read and annotation are on the same strand
	 *                          StrandStrategy.ANTISENSE retains all annotations where the read and annotation are on the opposite strand
	 *                          StrandStrategy.BOTH retains all annotations
	 * @param acceptedFunctions The set of accepted functions this gene can have to be retained
	 */
	public AbstractFunctionalDataProcessor(final StrandStrategy strandStrategy, final Collection<LocusFunction> acceptedFunctions) {
		this.strandStrategy = strandStrategy;
		this.acceptedFunctions = acceptedFunctions;
	}

	public AbstractFunctionalDataProcessor(final StrandStrategy strandStrategy, final LocusFunction[] acceptedFunctions) {
		this(strandStrategy, Arrays.asList(acceptedFunctions));
	}

	public StrandStrategy getStrandStrategy() {
		return strandStrategy;
	}

	public Collection<LocusFunction> getAcceptedFunctions() {
		return acceptedFunctions;
	}

	public abstract List<FunctionalData> getFilteredFunctionalData (final String [] genes, final String[] strands, final LocusFunction [] locusFunctions, final boolean readNegativeStrand);

	public abstract List<FunctionalData> filterToPreferredAnnotations(final Collection<FunctionalData> fdList);









}
