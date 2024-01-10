package org.broadinstitute.dropseqrna.annotation.functionaldata;

import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import picard.annotation.LocusFunction;

import java.util.Collection;
import java.util.List;

public class StarSoloFunctionalDataProcessor extends AbstractFunctionalDataProcessor implements FunctionalDataProcessorI {
    /**
     * Initialize a processor with a set of locus functions that are acceptable, and a strand strategy.
     *
     * @param strandStrategy    The strategy for filtering on read and gene strand:
     *                          StrandStrategy.SENSE retains all annotations where the read and annotation are on the same strand
     *                          StrandStrategy.ANTISENSE retains all annotations where the read and annotation are on the opposite strand
     *                          StrandStrategy.BOTH retains all annotations
     * @param acceptedFunctions The set of accepted functions this gene can have to be retained
     */
    public StarSoloFunctionalDataProcessor(StrandStrategy strandStrategy, Collection<LocusFunction> acceptedFunctions) {
        super(strandStrategy, acceptedFunctions);
    }

    public StarSoloFunctionalDataProcessor(StrandStrategy strandStrategy, LocusFunction[] acceptedFunctions) {
        super(strandStrategy, acceptedFunctions);
    }

    @Override
    public List<FunctionalData> getFilteredFunctionalData(String[] genes, String[] strands, LocusFunction[] locusFunctions, boolean readNegativeStrand) {
        return null;
    }

    @Override
    public List<FunctionalData> filterToPreferredAnnotations(Collection<FunctionalData> fdList) {
        return null;
    }
}
