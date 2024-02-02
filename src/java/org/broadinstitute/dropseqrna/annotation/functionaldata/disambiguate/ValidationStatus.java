package org.broadinstitute.dropseqrna.annotation.functionaldata.disambiguate;

import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalData;

public class ValidationStatus {

    private final GeneWithFunction star;


    private final GeneWithFunction dropseq;

    public ValidationStatus (GeneWithFunction star, GeneWithFunction dropseq) {
        this.dropseq=dropseq;
        this.star=star;
    }

    public boolean isValid () {
        // for antisense annotations, ignore the gene label.
        if (star.getCategory()== FunctionalData.Type.CODING_ANTISENSE && dropseq.getCategory()== FunctionalData.Type.CODING_ANTISENSE)
            return true;

        // otherwise, straight up comparison.
        boolean flag = star.equals(dropseq);
        return flag;

    }

    public GeneWithFunction getStar() {
        return star;
    }

    public GeneWithFunction getDropseq() {
        return dropseq;
    }



}
