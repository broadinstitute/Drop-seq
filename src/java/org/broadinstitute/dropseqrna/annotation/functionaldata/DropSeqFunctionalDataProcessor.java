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
public class DropSeqFunctionalDataProcessor extends AbstractFunctionalDataProcessor implements FunctionalDataProcessorI {


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
        super(strandStrategy, acceptedFunctions);
    }

    public DropSeqFunctionalDataProcessor(final StrandStrategy strandStrategy, final LocusFunction[] acceptedFunctions) {
        this(strandStrategy, Arrays.asList(acceptedFunctions));
    }

    /*
     * All lists have the same number of elements
     */
    private List<FunctionalData> convert (final String[] genes, final String[] strands, final LocusFunction[] locusFunctions) {
        List<FunctionalData> data = new ArrayList<FunctionalData>();
        if (getStrandStrategy() != null && genes.length != strands.length) {
            throw new IllegalArgumentException("Genes and strands must be of the same length.");
        }
        if (!getAcceptedFunctions().isEmpty() && genes.length != locusFunctions.length)
            throw new IllegalArgumentException("Genesand locus functions must be of the same length.");
        for (int i = 0; i < genes.length; i++) {
            FunctionalData sf = new FunctionalData(genes[i],
                    (getStrandStrategy() == null ? null : strands[i]),
                    (getAcceptedFunctions().isEmpty() ? null : locusFunctions[i]));
            data.add(sf);
        }
        return (data);
    }

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
        List<FunctionalData> result = convert(genes, strands, locusFunctions);
        if (getStrandStrategy() != null) {
            result = filterOnStrand(getStrandStrategy(), readNegativeStrand, result);
        }
        if (!getAcceptedFunctions().isEmpty()) {
            result = filterOnLocus(getAcceptedFunctions(), result);
        }
        // need to reduce result down to a unique set of genes/strands.
        // so if a gene is matched at UTR and CODING separately, that counts as 1 read.
        result = simplifyFD(result);
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

        return result;
    }

    private boolean isCodingUTR(final LocusFunction f) {
        return (f == LocusFunction.CODING || f == LocusFunction.UTR);
    }

    /**
     * Find the unique set of genes and set the best functional annotation for each gene.
     */
    private List<FunctionalData> simplifyFD(final List<FunctionalData> data) {
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
            FunctionalData result = new FunctionalData(l.get(0).getGene(), l.get(0).getStrand(), finalLF);
            resultList.add(result);
        }
        return (resultList);
    }

    /**
     * Return the set of functional annotations this object holds by the strand of a read and strand filtering strategy.
     * GENE_FUNCTION_STRAND.SENSE retains all annotations where the read and annotation are on the same strand
     * GENE_FUNCTION_STRAND.ANTISENSE retains all annotations where the read and annotation are on the opposite strand
     * GENE_FUNCTION_STRAND.BOTH retains all annotations
     *
     * @param strand             the filtering strategy for annotations
     * @param readNegativeStrand is the read on the negative strand?
     * @param originalData       the set of data to filter.
     * @return the subset of FunctionalData passing this filter.
     */
    private List<FunctionalData> filterOnStrand(final StrandStrategy strand, final boolean readNegativeStrand, final List<FunctionalData> originalData) {
        String readStrandString = Utils.strandToString(!readNegativeStrand);
        List<FunctionalData> result = new ArrayList<FunctionalData>();

        if (strand == StrandStrategy.SENSE)
            for (FunctionalData sf : originalData)
                if (sf.getStrand().equals(readStrandString))
                    result.add(sf);

        if (strand == StrandStrategy.ANTISENSE)
            for (FunctionalData sf : originalData)
                if (!sf.getStrand().equals(readStrandString))
                    result.add(sf);

        // NO OP
        if (strand == StrandStrategy.BOTH)
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
     * @param functions    A list of function types to retain.
     * @param originalData the set of data to filter.
     * @return the subset of FunctionalData passing this filter.
     */
    private List<FunctionalData> filterOnLocus(final Collection<LocusFunction> functions, final List<FunctionalData> originalData) {
        List<FunctionalData> result = new ArrayList<FunctionalData>();
        // collection data by gene.
        Map<String, List<FunctionalData>> originalDataByGene = new HashMap<String, List<FunctionalData>>();
        for (FunctionalData sf : originalData) {
            List<FunctionalData> fd = originalDataByGene.get(sf.getGene());
            if (fd == null) {
                fd = new ArrayList<FunctionalData>();
                originalDataByGene.put(sf.getGene(), fd);
            }
            fd.add(sf);
        }
        // test by gene
        for (String gene : originalDataByGene.keySet()) {
            List<FunctionalData> fd = originalDataByGene.get(gene);
            fd = filterOnLocusFunctionByGene(functions, fd);
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
        Set<LocusFunction> lfByGene = new HashSet<LocusFunction>();
        for (FunctionalData fd : originalDataByGene)
            lfByGene.add(fd.getLocusFunction());
        if (functions.containsAll(lfByGene))
            return originalDataByGene;
        return Collections.emptyList();

    }


}
