package org.broadinstitute.dropseqrna.annotation.functionaldata;

import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import picard.annotation.LocusFunction;

import java.util.Collection;
import java.util.List;

public interface FunctionalDataProcessorI {

    /**
     * Return the set of functional annotations this object holds by the strand of a read and strand filtering strategy,
     * and by the set of genes that overlap accepted locus functions.  If a gene has multiple functional annotations,
     * the "best" annotation is selected (see: AnnotationUtils.getLocusFunction()).
     *
     * @param genes An array of gene names
     * @param strands The strands of the genes
     * @param locusFunctions The locus function of each gene
     * @param readNegativeStrand Is the read these gene annotations come from on the negative strand?
     * @return A list of FunctionalData that passes criteria.
     */
    public List<FunctionalData> getFilteredFunctionalData (final String [] genes, final String[] strands, final LocusFunction[] locusFunctions, final boolean readNegativeStrand);

    /**
     * Across all the functional data, find the annotation type with the highest priority.  This priority order is
     * implementation specific for a particular strategy.
     * Restrict functional annotation data to that subset.
     * @param fdList A list of functional data to be filtered
     * @return A subset of the input list of functional data.
     */
    public List<FunctionalData> filterToPreferredAnnotations(final Collection<FunctionalData> fdList);

    public StrandStrategy getStrandStrategy();

    public Collection<LocusFunction> getAcceptedFunctions();
}
