package org.broadinstitute.dropseqrna.eqtl;

public enum DonorMergeStrategy {	
        Sum("For each gene, sum the UMIs of each cell into a single MetaCell"),
        Mean("For each gene, calculate the mean number of UMIs across cells into a single MetaCell"),
        Median("For each gene, calculate the median number of UMIs across cells into a single MetaCell");

        public final String description;

        DonorMergeStrategy(final String description) {
            this.description = description;
        }

        /** Gets the description of the strategy. */
        public String getDescription() {
            return description;
        }    
}
