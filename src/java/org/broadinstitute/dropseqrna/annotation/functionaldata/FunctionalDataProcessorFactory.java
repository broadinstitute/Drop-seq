package org.broadinstitute.dropseqrna.annotation.functionaldata;

import org.apache.commons.lang3.NotImplementedException;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import picard.annotation.LocusFunction;

import java.util.Collection;

public class FunctionalDataProcessorFactory {
    public static FunctionalDataProcessorI getFunctionalDataProcessor(final StrandStrategy strandStrategy, final Collection<LocusFunction> acceptedFunctions, FunctionalDataProcessorStrategy strategy) {
        if (strategy== FunctionalDataProcessorStrategy.DROPSEQ)
            return new DropSeqFunctionalDataProcessor(strandStrategy, acceptedFunctions);
        if (strategy== FunctionalDataProcessorStrategy.STARSOLO)
            return new StarSoloFunctionalDataProcessor(strandStrategy, acceptedFunctions);
        throw new NotImplementedException("Strategy" + strategy +"not implemented!");
    }

}
