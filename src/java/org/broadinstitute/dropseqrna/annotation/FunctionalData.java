package org.broadinstitute.dropseqrna.annotation;

import picard.annotation.LocusFunction;

/**
 * A simple class to hold the tuple of information about a read's gene, strand, function type for a single gene.
 * @author nemesh
 *
 */
public class FunctionalData {

		private String gene;
		private String strand;
		private LocusFunction locusFunction;

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
}
