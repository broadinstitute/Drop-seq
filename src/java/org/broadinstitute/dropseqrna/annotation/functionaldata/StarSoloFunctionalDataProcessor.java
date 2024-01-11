package org.broadinstitute.dropseqrna.annotation.functionaldata;

import org.broadinstitute.dropseqrna.annotation.AnnotationUtils;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import picard.annotation.LocusFunction;

import java.util.*;

public class StarSoloFunctionalDataProcessor implements FunctionalDataProcessorI {

    private final DataProcessorUtils util;
    private final PriorityScoreI priority = new StarSoloPriorityScore();

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
        util=new DataProcessorUtils(strandStrategy, acceptedFunctions);
    }

    public StarSoloFunctionalDataProcessor(StrandStrategy strandStrategy, LocusFunction[] acceptedFunctions) {
        this(strandStrategy, Arrays.asList(acceptedFunctions));
    }




    /**
     * Return the set of functional annotations this object holds by the strand of a read and strand filtering strategy,
     * and by the set of genes that overlap accepted locus functions.  If a gene has multiple functional annotations,
     * the "best" annotation is selected (see: AnnotationUtils.getLocusFunction()).
     * The gene functions are not filtered by the strand of the read to preserve "antisense" coding annotations.
     */
    public List<FunctionalData> getFilteredFunctionalData(String[] genes, String[] strands, LocusFunction[] locusFunctions, boolean readNegativeStrand) {
        List<FunctionalData> result = FunctionalData.buildFD(genes, strands, locusFunctions, this.util.getStrandStrategy(), this.util.getAcceptedFunctions(),readNegativeStrand);
        if (!getAcceptedFunctions().isEmpty()) {
            result = util.filterOnLocus(result);
        }
        // need to reduce result down to a unique set of genes/strands.
        // so if a gene is matched at UTR and CODING separately, that counts as 1 read.
        result = util.simplifyFD(result);
        return (result);
    }


    /**
     * Filter from multiple gene/functional annotations to a single annotation using the STARsolo priority.
     * In this pattern, the off-strand functional annotation has been retained by getFilteredFunctionalData,
     * so off-strand exonic/UTR regions can compete (and outscore) on-strand introic regions.
     *
     * @param fdList A list of functional data to be filtered
     * @return
     */
    public List<FunctionalData> filterToPreferredAnnotations(Collection<FunctionalData> fdList) {
        List<FunctionalData> result = util.filterToPreferredAnnotations(fdList, priority);
        return result;
    }

    public StrandStrategy getStrandStrategy() {
        return this.util.getStrandStrategy();
    }

    @Override
    public Collection<LocusFunction> getAcceptedFunctions() {
        return this.util.getAcceptedFunctions();
    }

}
