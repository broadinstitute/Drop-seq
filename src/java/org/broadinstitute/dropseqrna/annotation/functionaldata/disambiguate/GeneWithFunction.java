package org.broadinstitute.dropseqrna.annotation.functionaldata.disambiguate;

import com.google.common.base.Objects;
import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalData;

public class GeneWithFunction {

    private final String gene;
    private final FunctionalData.Type category;

    public GeneWithFunction (String gene, FunctionalData.Type category) {
        this.gene=gene;
        this.category=category;
    }

    public String getGene() {
        return gene;
    }

    public FunctionalData.Type getCategory() {
        return category;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        GeneWithFunction that = (GeneWithFunction) o;
        return Objects.equal(gene, that.gene) && category == that.category;
    }

    @Override
    public int hashCode() {
        return Objects.hashCode(gene, category);
    }


    @Override
    public String toString() {
        return gene + " [" + category+"]";
    }
}
