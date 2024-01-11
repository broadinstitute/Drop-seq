package org.broadinstitute.dropseqrna.annotation.functionaldata;

import org.broadinstitute.dropseqrna.annotation.AnnotationUtils;
import org.broadinstitute.dropseqrna.barnyard.Utils;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import picard.annotation.LocusFunction;

import java.util.*;

/**
 * Given lists of gene names, strands of the gene, and functional annotations of the genes,
 * decide which genes/strands/functions should be retained.
 *
 * @author nemesh
 */
public class DropSeqFunctionalDataProcessor implements FunctionalDataProcessorI {

    private final DataProcessorUtils util;
    private final PriorityScoreI priority = new DropSeqPriorityScore();
    /**
     * Initialize a processor with a set of locus functions that are acceptable, and a strand strategy.
     *
     * @param strandStrategy    The strategy for filtering on read and gene strand:
     *                          StrandStrategy.SENSE retains all annotations where the read and annotation are on the same strand
     *                          StrandStrategy.ANTISENSE retains all annotations where the read and annotation are on the opposite strand
     *                          StrandStrategy.BOTH retains all annotations
     * @param acceptedFunctions The set of accepted functions this gene can have to be retained
     */
    public DropSeqFunctionalDataProcessor(final StrandStrategy strandStrategy, final Collection<LocusFunction> acceptedFunctions) {
        util=new DataProcessorUtils(strandStrategy, acceptedFunctions);
    }

    public DropSeqFunctionalDataProcessor(final StrandStrategy strandStrategy, final LocusFunction[] acceptedFunctions) {
        this(strandStrategy, Arrays.asList(acceptedFunctions));
    }

    /*
     * All lists have the same number of elements
     */

    /**
     * Return the set of functional annotations this object holds by the strand of a read and strand filtering strategy,
     * and by the set of genes that overlap accepted locus functions.  If a gene has multiple functional annotations,
     * the "best" annotation is selected (see: AnnotationUtils.getLocusFunction()).
     *
     * @param readNegativeStrand Is the read on the negative strand?
     * @return A list of FunctionalData that passes criteria.
     */
    @Override
	public List<FunctionalData> getFilteredFunctionalData(final String[] genes, final String[] strands, final LocusFunction[] locusFunctions, final boolean readNegativeStrand) {
        List<FunctionalData> result = FunctionalData.buildFD(genes, strands, locusFunctions, this.util.getStrandStrategy(), this.util.getAcceptedFunctions(), readNegativeStrand);
        if (util.getStrandStrategy() != null) {
            result = util.filterOnStrand(readNegativeStrand, result);
        }
        if (!getAcceptedFunctions().isEmpty()) {
            result = util.filterOnLocus(result);
        }
        // need to reduce result down to a unique set of genes/strands.
        // so if a gene is matched at UTR and CODING separately, that counts as 1 read.
        result = util.simplifyFD(result);
        return (result);
    }


    /**
     * Across all the functional data, find the annotation type with the highest priority (coding/utr>intron>intergenic)
     * Restrict functional annotation data to that subset.
     * The desire here is to take reads that are assigned to both the coding region of gene A and the intronic region of gene B
     * and filter out gene B, resulting in a read that is unambiguously assigned.
     *
     * @param fdList A list of functional data to be filtered
     * @return A subset of the input list of functional data.
     */
    @Override
	public List<FunctionalData> filterToPreferredAnnotations(final Collection<FunctionalData> fdList) {
        /*
        boolean hasCoding = false;

        for (FunctionalData fd : fdList) {
            if (isCodingUTR(fd.getLocusFunction())) {
                hasCoding = true;
                break;
            }
        }

        List<FunctionalData> result = new ArrayList<FunctionalData>();
        for (FunctionalData fd : fdList) {
            boolean readIsCoding = isCodingUTR(fd.getLocusFunction());
            if (!hasCoding || (readIsCoding)) {
                result.add(fd);
            }
        }
        */
        List<FunctionalData> result = util.filterToPreferredAnnotations(fdList, priority);
        return result;
    }

    @Override
    public StrandStrategy getStrandStrategy() {
        return this.util.getStrandStrategy();
    }

    @Override
    public Collection<LocusFunction> getAcceptedFunctions() {
        return this.util.getAcceptedFunctions();
    }

    private boolean isCodingUTR(final LocusFunction f) {
        return (f == LocusFunction.CODING || f == LocusFunction.UTR);
    }






}
