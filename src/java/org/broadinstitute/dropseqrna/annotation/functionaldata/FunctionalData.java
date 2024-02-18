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

	private final boolean isAmbiguous;
	private Type type;

	public FunctionalData (final String gene, final String geneStrand, final LocusFunction locusFunction, final String readStrand) {
		this.gene=gene;
		this.geneStrand =geneStrand;
		this.locusFunction=locusFunction;
		this.readStrand=readStrand;
		this.isAmbiguous=false;
	}

	public FunctionalData (final String gene, final String geneStrand, final LocusFunction locusFunction, final String readStrand, boolean isAmbiguous) {
		this.gene=gene;
		this.geneStrand =geneStrand;
		this.locusFunction=locusFunction;
		this.readStrand=readStrand;
		this.isAmbiguous=isAmbiguous;
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
					Utils.negativeStrandToString(readNegativeStrand));
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

	/**
	 * Test if the read strand and gene strand match
	 * @return True if the read strand == gene strand
	 */
	public boolean isSense () {
		return (this.getGeneStrand().equals(this.readStrand));
	}

	public LocusFunction getLocusFunction () {
		return this.locusFunction;
	}

	@Override
	public String toString () {
		return String.format("%s %s %s %s %s" ,
				java.util.Objects.toString(gene, "null"),
				java.util.Objects.toString(geneStrand, "null"),
				java.util.Objects.toString(locusFunction, "null"),
				java.util.Objects.toString(readStrand, "null"),
				java.util.Objects.toString(getType(), "null")
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

	/**
	 * Test if two FunctionalData objects have the same interpretation - the same gene and type.
	 * For example, for a metagene, the one copy might align on the + strand, and one on the negative strand
	 * They both interpret to CODING_SENSE, but equals will say they are different.
	 * @param other
	 * @return
	 */
	public boolean sameGeneAndType (FunctionalData other) {
		return this.getType()==other.getType() & this.getGene().equals(other.getGene());
	}

	/**
	 * A convenience method to get the type this functional data encodes.
	 * This combines the gene strand, read strand, and locus function to yield a single classification
	 * This simplifies CODING and UTR into CODING.
	 * @return
	 */
	public Type getType() {
		if (this.type == null) {
			this.type = calculateType();
		}
		return this.type;
	}

	private Type calculateType() {
		if (this.locusFunction == LocusFunction.INTERGENIC && this.isAmbiguous)
			return Type.AMBIGUOUS;

		if (this.locusFunction == LocusFunction.INTERGENIC) {
			return Type.INTERGENIC;
		}

		if (this.isSense()) {
			switch (this.locusFunction) {
				case CODING:
				case UTR:
					return Type.CODING_SENSE;
				case INTRONIC:
					return Type.INTRONIC_SENSE;
			}
		} else {
			switch (this.locusFunction) {
				case CODING:
				case UTR:
					return Type.CODING_ANTISENSE;
				case INTRONIC:
					return Type.INTRONIC_ANTISENSE;
			}
		}

		// Anything else
		return Type.INTERGENIC;
	}

	public enum Type {
		CODING_SENSE,
		CODING_ANTISENSE,
		INTRONIC_SENSE,
		INTRONIC_ANTISENSE,
		INTERGENIC,
		AMBIGUOUS
	}
}
