package org.broadinstitute.dropseqrna.annotation.functionaldata;

import org.broadinstitute.dropseqrna.annotation.AnnotationUtils;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import picard.annotation.LocusFunction;

import java.util.*;

public class DataProcessorUtils {

    private final StrandStrategy strandStrategy;
    private final Collection<LocusFunction> acceptedFunctions;

    /**
     * Initialize a processor with a set of locus functions that are acceptable, and a strand strategy.
     *
     * @param strandStrategy    The strategy for filtering on read and gene strand:
     *                          StrandStrategy.SENSE retains all annotations where the read and annotation are on the same strand
     *                          StrandStrategy.ANTISENSE retains all annotations where the read and annotation are on the opposite strand
     *                          StrandStrategy.BOTH retains all annotations
     * @param acceptedFunctions The set of accepted functions this gene can have to be retained
     */
    public DataProcessorUtils (final StrandStrategy strandStrategy, final Collection<LocusFunction> acceptedFunctions) {
        this.strandStrategy = strandStrategy;
        this.acceptedFunctions = acceptedFunctions;
    }

    public StrandStrategy getStrandStrategy() {
        return this.strandStrategy;
    }

    public Collection<LocusFunction> getAcceptedFunctions() {
        return this.acceptedFunctions;
    }

    /**
     * Return the set of functional annotations this object holds by the strand of a read and strand filtering strategy.
     * GENE_FUNCTION_STRAND.SENSE retains all annotations where the read and annotation are on the same strand
     * GENE_FUNCTION_STRAND.ANTISENSE retains all annotations where the read and annotation are on the opposite strand
     * GENE_FUNCTION_STRAND.BOTH retains all annotations
     *
     * @param readNegativeStrand is the read on the negative strand?
     * @param originalData       the set of data to filter.
     * @return the subset of FunctionalData passing this filter.
     */
    public List<FunctionalData> filterOnStrand(final boolean readNegativeStrand, final List<FunctionalData> originalData) {

        List<FunctionalData> result = new ArrayList<>();

        if (this.strandStrategy == StrandStrategy.SENSE)
            for (FunctionalData sf : originalData)
                if (sf.isGeneNegativeStrand()==readNegativeStrand)
                    result.add(sf);

        if (this.strandStrategy == StrandStrategy.ANTISENSE)
            for (FunctionalData sf : originalData)
                if (sf.isGeneNegativeStrand()!=readNegativeStrand)
                    result.add(sf);

        // NO OP
        if (this.strandStrategy == StrandStrategy.BOTH)
            return (originalData);

        return result;
    }

    /**
     * Return the set of functional annotations this object holds by the list of functional types that are acceptable.
     * <p>
     * For a single gene's functional data, are all of the locus functions in the collection of accepted locus function?
     * IE: if the locus functions accepted are UTR/CODING and the gene has CODING or UTR or CODING+UTR, accept these functions.
     * If the gene has CODING+UTR+INTRON, reject the entire gene.  In this case, return an empty list.
     *
     * @param originalData the set of data to filter.
     * @return the subset of FunctionalData passing this filter.
     */
    public List<FunctionalData> filterOnLocus(final List<FunctionalData> originalData) {
        List<FunctionalData> result = new ArrayList<>();
        // collection data by gene.
        Map<String, List<FunctionalData>> originalDataByGene = new HashMap<>();
        for (FunctionalData sf : originalData) {
            List<FunctionalData> fd = originalDataByGene.computeIfAbsent(sf.getGene(), k -> new ArrayList<>());
            fd.add(sf);
        }
        // test by gene
        for (String gene : originalDataByGene.keySet()) {
            List<FunctionalData> fd = originalDataByGene.get(gene);
            fd = filterOnLocusFunctionByGene(acceptedFunctions, fd);
            result.addAll(fd);
        }
        return result;
    }

    /**
     * For a single gene's functional data, are all of the locus functions in the collection of accepted locus function?
     * IE: if the locus functions accepted are UTR/CODING and the gene has CODING or UTR or CODING+UTR, accept these functions.
     * If the gene has CODING+UTR+INTRON, reject the entire gene.  In this case, return an empty list.
     *
     * @param functions          The list of accepted locus functions.
     * @param originalDataByGene The functional data for a single gene.
     * @return filtered list of FunctionalData
     */
    private List<FunctionalData> filterOnLocusFunctionByGene(final Collection<LocusFunction> functions, final List<FunctionalData> originalDataByGene) {
        Set<LocusFunction> lfByGene = new HashSet<>();
        for (FunctionalData fd : originalDataByGene)
            lfByGene.add(fd.getLocusFunction());
        if (functions.containsAll(lfByGene))
            return originalDataByGene;
        return Collections.emptyList();

    }

    /**
     * Find the unique set of genes and set the best functional annotation for each gene.
     */
    public List<FunctionalData> simplifyFD(final List<FunctionalData> data) {
        // map each gene to the list of functional data, then reduce each gene.
        Map<String, List<FunctionalData>> map = new HashMap<String, List<FunctionalData>>();
        for (FunctionalData fd : data) {
            List<FunctionalData> l = map.get(fd.getGene());
            if (l == null) {
                l = new ArrayList<FunctionalData>();
            }
            l.add(fd);
            map.put(fd.getGene(), l);
        }

        List<FunctionalData> resultList = new ArrayList<FunctionalData>();

        // process each list.
        for (String gene : map.keySet()) {
            List<FunctionalData> l = map.get(gene);
            final LocusFunction finalLF;
            if (getAcceptedFunctions().isEmpty()) {
                finalLF = null;
            } else {
                List<LocusFunction> lf = new ArrayList<LocusFunction>();
                for (FunctionalData fd : l)
                    lf.add(fd.getLocusFunction());
                finalLF = AnnotationUtils.getInstance().getLocusFunction(lf, false);
            }
            // any instance of list l is fine.
            FunctionalData result = new FunctionalData(l.getFirst().getGene(), l.getFirst().isGeneNegativeStrand(), finalLF, l.getFirst().isReadNegativeStrand());
            resultList.add(result);
        }
        return (resultList);
    }

    /***
     * Filter a list of functional data to the set that have the highest priority, as defined by the PriorityScore object.
     * This is useful to filter out lower level functional data, such as when an intronic and coding gene are both captured
     * by the same read - the gene(s) encoding the coding region are kept, and the lower priority introns removed.
     *
     * @param fdList A list of functional data to filter
     * @param priorityScore An method to score functional data objects to promote some as more important.
     * @return
     */
    public List<FunctionalData> filterToPreferredAnnotations(final Collection<FunctionalData> fdList, PriorityScoreI priorityScore) {
        int bestScore = Integer.MAX_VALUE;
        List<FunctionalData> result = new ArrayList<>();

        for (FunctionalData fd: fdList) {
            int score = priorityScore.getScore(fd);
            // if the score is better, clear the list and start remembering functional data
            if (score < bestScore) {
                bestScore = score;
                result.clear();
            }
            if (score==bestScore) {
                result.add(fd);
            }
        }
        return result;

    }

}



