package org.broadinstitute.dropseqrna.annotation.functionaldata.disambiguate;

import org.broadinstitute.dropseqrna.utils.readiterators.GeneFunctionProcessor;

public class ClassifyFunctionalDataBase {
    final String missingGeneLabel;
    final GeneFunctionProcessor gfp;

    public static final String DEFAULT_MISSING_GENE_LABEL="-";

    public ClassifyFunctionalDataBase(GeneFunctionProcessor gfp, String missingGeneLabel) {
        this.gfp=gfp;
        this.missingGeneLabel=missingGeneLabel;
    }

    public ClassifyFunctionalDataBase(GeneFunctionProcessor gfp) {
        this(gfp, DEFAULT_MISSING_GENE_LABEL);
    }

}
