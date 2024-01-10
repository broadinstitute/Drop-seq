package org.broadinstitute.dropseqrna.annotation.functionaldata;

import com.google.common.base.Objects;
import picard.annotation.LocusFunction;

/**
 * A simple class to hold the tuple of information about a read's gene, strand, function type for a single gene.
 * @author nemesh
 *
 */
public class FunctionalData {

		private final String gene;
		private final String strand;
		private final LocusFunction locusFunction;

		public FunctionalData (final String gene, final String strand, final LocusFunction locusFunction) {
			this.gene=gene;
			this.strand=strand;
			this.locusFunction=locusFunction;
		}

		public String getStrand () {
			return this.strand;
		}

		public String getGene () {
			return this.gene;
		}

		public LocusFunction getLocusFunction () {
			return this.locusFunction;
		}

		@Override
		public String toString () {
			return (this.gene + " " + this.strand + " " + this.locusFunction.name());
		}

	@Override
	public boolean equals(Object o) {
		if (this == o) return true;
		if (o == null || getClass() != o.getClass()) return false;
		FunctionalData that = (FunctionalData) o;
		return Objects.equal(gene, that.gene) && Objects.equal(strand, that.strand) && locusFunction == that.locusFunction;
	}

	@Override
	public int hashCode() {
		return Objects.hashCode(gene, strand, locusFunction);
	}
}
