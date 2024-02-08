package org.broadinstitute.dropseqrna.annotation.functionaldata.disambiguate;

import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalData;
import org.broadinstitute.dropseqrna.utils.readiterators.GeneFunctionProcessor;

import java.util.List;

public class ClassifyStarSoloFunctionalData extends ClassifyFunctionalDataBase {

    public ClassifyStarSoloFunctionalData(GeneFunctionProcessor gfp, String missingGeneLabel) {
        super(gfp, missingGeneLabel);
    }

    public ClassifyStarSoloFunctionalData(GeneFunctionProcessor gfp) {
        super(gfp);
    }

    public GeneWithFunction convert (FunctionalData fd) {
        GeneWithFunction result = new GeneWithFunction(fd.getGene(), fd.getType());
        return result;
    }


}
