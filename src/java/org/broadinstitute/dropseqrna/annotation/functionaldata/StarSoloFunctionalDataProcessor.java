package org.broadinstitute.dropseqrna.annotation.functionaldata;

import org.broadinstitute.dropseqrna.annotation.AnnotationUtils;
import org.broadinstitute.dropseqrna.barnyard.Utils;
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


    public List<FunctionalData> getFilteredFunctionalData(List <FunctionalData> fdList) {
        if (!getAcceptedFunctions().isEmpty()) {
            fdList = util.filterOnLocus(fdList);
        }
        // need to reduce result down to a unique set of genes/strands.
        // so if a gene is matched at UTR and CODING separately, that counts as 1 read.
        fdList = util.simplifyFD(fdList);
        return (fdList);
    }

    /**
     * Filter from multiple gene/functional annotations to a single annotation using the STARsolo priority.
     * In this pattern, the off-strand functional annotation has been retained by getFilteredFunctionalData,
     * so off-strand exonic/UTR regions can compete (and outscore) on-strand introic regions.
     *
     * @param fdList A list of functional data to be filtered
     * @return
     */
    public List<FunctionalData> filterToPreferredAnnotations(List<FunctionalData> fdList) {
        List<FunctionalData> result = util.filterToPreferredAnnotations(fdList, priority);
        result = util.filterOnStrand(result);
        return result;
    }

    public StrandStrategy getStrandStrategy() {
        return this.util.getStrandStrategy();
    }

    @Override
    public Collection<LocusFunction> getAcceptedFunctions() {
        return this.util.getAcceptedFunctions();
    }

    @Override
    public PriorityScoreI getPriority() {
        return this.priority;
    }

    /**
     * Convert the Starsolo read annotations to a FunctionalData object.
     * @param geneName The name of the gene
     * @param priority The gene function priority as interpreted by StarSolo.  Ranges from 1 to 8, with
     *                     higher numbers having higher priority.
     * @param numGenes How many genes does this read map to at the priority score.  If more than one,
     *                 the read is considered ambiguous, which is considered as Intergenic by Starsolo to
     *                 encode "not counted".  We encode this as AMBIGUOUS instead of INTERGENIC to differentiate
     *                 between a read that maps to an intergenic region and a read that maps to multiple genes.
     * @param readNegativeStrand True if the read is on the negative strand
     * @return A FunctionalData object encoding the results of this read.
     */
    public static FunctionalData getFunctionalData (String geneName, int priority, int numGenes, boolean readNegativeStrand) {
        LocusFunction lf = convertPriority(priority, numGenes);
        String readStrand = Utils.negativeStrandToString(readNegativeStrand);

        // we have to infer the strand of the gene from the read and the priority.
        // read matches gene for priority 1,3,5.  Strand opposite read 2,4,6.
        String geneStrand = readStrand;
        if (priority % 2==0)
            geneStrand = Utils.negativeStrandToString(!readNegativeStrand);
        boolean isAmbiguous = numGenes >1;
        FunctionalData fd = new FunctionalData(geneName, geneStrand, lf, readStrand, isAmbiguous);
        return (fd);
    }

    /**
     * Convert the StarSolo numeric priority score to a LocusFunction
     * This handles the type of locus function, but not the interaction of strand of gene and read.
     * @param priority The priority score to convert
     * @return The LocusFunction encoded.
     */
    public static LocusFunction convertPriority (int priority, int numGenes) {
        if (numGenes> 1)
            return LocusFunction.INTERGENIC;
        LocusFunction result = switch (priority) {
            case 1,2,3,4-> {
                yield LocusFunction.CODING;
            }
            case 5,6 -> {
                yield LocusFunction.INTRONIC;
            }
            default ->LocusFunction.INTERGENIC;
        };
        return result;
    }

}
