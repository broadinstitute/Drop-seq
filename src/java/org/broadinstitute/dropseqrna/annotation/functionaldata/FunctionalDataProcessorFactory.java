package org.broadinstitute.dropseqrna.annotation.functionaldata;

import org.apache.commons.lang3.NotImplementedException;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import picard.annotation.LocusFunction;

import java.util.Collection;

public class FunctionalDataProcessorFactory {
    public static FunctionalDataProcessorI getFunctionalDataProcessor(final StrandStrategy strandStrategy, final Collection<LocusFunction> acceptedFunctions, FunctionalDataProcessorStrategyEnum strategy) {
        if (strategy==FunctionalDataProcessorStrategyEnum.DROPSEQ)
            return new DropSeqFunctionalDataProcessor(strandStrategy, acceptedFunctions);
        if (strategy==FunctionalDataProcessorStrategyEnum.STARSOLO)
            return new StarSoloFunctionalDataProcessor(strandStrategy, acceptedFunctions);
        throw new NotImplementedException("Strategy" + strategy +"not implemented!");
    }

}
