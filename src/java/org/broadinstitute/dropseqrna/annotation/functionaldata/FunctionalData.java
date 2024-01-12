package org.broadinstitute.dropseqrna.annotation.functionaldata;

import com.google.common.base.Objects;
import org.broadinstitute.dropseqrna.barnyard.Utils;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import picard.annotation.LocusFunction;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * A simple class to hold the tuple of information about a read's gene, strand, function type for a single gene.
 * @author nemesh
 *
 */
public class FunctionalData {

	private final String gene;
	private final String geneStrand;
	private final LocusFunction locusFunction;

	private final String readStrand;

	public FunctionalData (final String gene, final String geneStrand, final LocusFunction locusFunction, final String readStrand) {
		this.gene=gene;
		this.geneStrand =geneStrand;
		this.locusFunction=locusFunction;
		this.readStrand=readStrand;
	}



	/**
	 * Construct a list of functional data objects from a set of raw inputs.
	 * This validates that the genes/strands/locusfunctions are of the same length
	 * @param genes An array of gene names
	 * @param strands An array of strands for the genes
	 * @param locusFunctions An array of locus functions for the genes
	 * @param strandStrategy The strategy used to interpret the strand.  If null, the strand value is set to null
	 * @param acceptedFunctions The accepted functions for the analysis.  If empty, the function for each FD is set to null.
	 * @return A list of FunctionData objects.
	 */
	public static List<FunctionalData> buildFD (final String[] genes, final String[] strands, final LocusFunction[] locusFunctions,
										 final StrandStrategy strandStrategy, final Collection<LocusFunction> acceptedFunctions,
												final boolean readNegativeStrand) {
		List<FunctionalData> data = new ArrayList<FunctionalData>();
		if (strandStrategy != null && genes.length != strands.length) {
			throw new IllegalArgumentException("Genes and strands must be of the same length.");
		}
		if (!acceptedFunctions.isEmpty() && genes.length != locusFunctions.length)
			throw new IllegalArgumentException("Genes and locus functions must be of the same length.");
		for (int i = 0; i < genes.length; i++) {
			FunctionalData sf = new FunctionalData(genes[i],
					(strandStrategy == null ? null : strands[i]),
					(acceptedFunctions.isEmpty() ? null : locusFunctions[i]),
					Utils.strandToString(readNegativeStrand));
			data.add(sf);
		}
		return (data);
	}

	public String getGeneStrand() {
		return this.geneStrand;
	}

	public String getReadStrand() {
		return readStrand;
	}

	public String getGene () {
		return this.gene;
	}

	public LocusFunction getLocusFunction () {
		return this.locusFunction;
	}

	@Override
	public String toString () {
		return String.format("%s %s %s",
				java.util.Objects.toString(gene, "null"),
				java.util.Objects.toString(geneStrand, "null"),
				java.util.Objects.toString(locusFunction, "null"),
				java.util.Objects.toString(readStrand, "null")
				);
	}

	@Override
	public boolean equals(Object o) {
		if (this == o) return true;
		if (o == null || getClass() != o.getClass()) return false;
		FunctionalData that = (FunctionalData) o;
		return Objects.equal(gene, that.gene) && Objects.equal(geneStrand, that.geneStrand) && locusFunction == that.locusFunction;
	}

	@Override
	public int hashCode() {
		return Objects.hashCode(gene, geneStrand, locusFunction);
	}
}
