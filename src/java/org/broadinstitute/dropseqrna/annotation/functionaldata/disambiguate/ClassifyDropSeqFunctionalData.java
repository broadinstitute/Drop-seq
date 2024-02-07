package org.broadinstitute.dropseqrna.annotation.functionaldata.disambiguate;

import com.google.api.gax.rpc.UnimplementedException;
import org.broadinstitute.dropseqrna.annotation.functionaldata.DataProcessorUtils;
import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalData;
import org.broadinstitute.dropseqrna.utils.readiterators.GeneFunctionProcessor;
import picard.annotation.LocusFunction;

import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;


/**
 * Simplify a list of functional data for a read to a simple structure containing
 * a gene and a function.
 *
 * This deals with lists of raw functional annotation data where there may be multiple
 * overlapping gene annotations, the result may be ambiguous, etc.
 */
public class ClassifyDropSeqFunctionalData extends ClassifyFunctionalDataBase{

    private DataProcessorUtils util;
    public ClassifyDropSeqFunctionalData(GeneFunctionProcessor gfp, String missingGeneLabel) {
        super(gfp, missingGeneLabel);
        util = new DataProcessorUtils(gfp.getFdp().getStrandStrategy(), gfp.getFdp().getAcceptedFunctions());
    }

    public ClassifyDropSeqFunctionalData(GeneFunctionProcessor gfp) {
        this(gfp, DEFAULT_MISSING_GENE_LABEL);
    }

    /**
     * Return the finalized gene and function for this list of dropseq functional annotations
     *
     * This differs slightly from the STAR implementation, as we retain the gene name for the antisense encoding.
     * @param fdList The input raw list of functional annotations for this read
     * @return The interpreted gene and function.
     */
    public GeneWithFunction convert (List<FunctionalData> fdList) {
        // begin testing conditions.
        // as each condition is tested, if a classification is reached return it.

        // to reduce headaches between CODING and UTR (, simplify all UTR to CODING
        List<FunctionalData> fdList2 = fdList.stream()
                .map(this::simplifyCoding)
                .collect(Collectors.toList());
        fdList=fdList2;

        // INTERGENIC - dropseq will contain an intergenic portion of the read, Star will flag the read as intergenic.
        boolean isIntergenicDropSeq=isIntergenic(fdList);
        if (isIntergenicDropSeq)
            return new GeneWithFunction(missingGeneLabel, FunctionalData.Type.INTERGENIC);

        // AMBIGUOUS - multiple genes assigned, no gene name.
        boolean isAmbiguous = isDropSeqAmbiguous(fdList);
        if (isAmbiguous)
            return new GeneWithFunction(missingGeneLabel, FunctionalData.Type.AMBIGUOUS);

        // ANTISENSE CODING
        FunctionalData antiSenseCoding = getDropSeqAntisenseCoding(fdList);
        if (antiSenseCoding!=null)
            return new GeneWithFunction(antiSenseCoding.getGene(), FunctionalData.Type.CODING_ANTISENSE);

        // ANTISENSE INTRONIC
        FunctionalData antiSenseIntronic = getDropSeqAntisenseIntronic(fdList);
        if (antiSenseIntronic!=null)
            return new GeneWithFunction(antiSenseIntronic.getGene(), FunctionalData.Type.INTRONIC_ANTISENSE);

        // SENSE (coding+intronic)
        FunctionalData sense = getDropSeqSense(fdList);
        if (sense!=null)
            return new GeneWithFunction(sense.getGene(), sense.getType());

        // sigh
        throw new UnsupportedOperationException("Data not classified.  You should not get here.");
    }



    boolean isIntergenic (List<FunctionalData> fdList) {
        List<FunctionalData> filtered = filterFunctionalData (fdList);
        boolean isIntergenicDropSeq=containsFunctionalAnnotation(fdList, LocusFunction.INTERGENIC);
        return filtered.isEmpty() && isIntergenicDropSeq;
    }

    /**
     * Ambiguous assignments have multiple functional data annotations to different genes at the same highest priority.
     * This counts as ambiguous reads that might not be counted (for example multiple genes associated to intronic antisense)
     * @param fdList
     * @return
     */
    boolean isDropSeqAmbiguous (List<FunctionalData> fdList) {
        // first test if there's one or more functions that pass.
        List<FunctionalData> filtered = filterFunctionalData (fdList);
        int countUniqueGenes = filtered.stream().map(FunctionalData::getGene).collect(Collectors.toSet()).size();
        // if there are multiple on-strand genes, it's ambiguous.
        if (countUniqueGenes>=2)
            return true;
        // if there's a single on-strand gene, it's not ambiguous
        if (countUniqueGenes==1)
            return false;
        // if there are multiple off-strand genes.  This occurs when the number of filtered genes is 0.
        filtered = filterFunctionalDataWithoutStrand (fdList);
        countUniqueGenes = filtered.stream().map(FunctionalData::getGene).collect(Collectors.toSet()).size();
        return countUniqueGenes > 1;
    }

    /**
     * Antisense has a filtered size of 0 and at least one fd is antisense.
     * The size of 0 indicates that there is no better coding annotation.
     * @param fdList An unfiltered list of functional data
     * @return True if there filtered data has no gene, and there is at least one antisense annotation.
     */
    FunctionalData getDropSeqAntisenseCoding (List<FunctionalData> fdList) {
        return (getDropSeqAntisense(fdList, FunctionalData.Type.CODING_ANTISENSE));
    }

    FunctionalData getDropSeqAntisenseIntronic (List<FunctionalData> fdList) {
        return (getDropSeqAntisense(fdList, FunctionalData.Type.INTRONIC_ANTISENSE));
    }

    private FunctionalData getDropSeqAntisense (List<FunctionalData> fdList, FunctionalData.Type type) {
        List<FunctionalData> filtered = filterFunctionalData (fdList);

        // simplify before testing antisense.
        fdList = util.simplifyFD(fdList);
        List<FunctionalData> antisense = fdList.stream().filter(x -> x.getType()==type).toList();
        // simplify multiple annotations on the same gene.

        if (filtered.isEmpty() && antisense.size()==1) {
            return antisense.getFirst();
        }
        return null;
    }

    /**
     * The unambiguous intrepretation of a list of functional annotations.
     * Assumes that other tests have been run first (ambiguous, antisense, intergenic)
     * @param fdList
     * @return
     */
    FunctionalData getDropSeqSense (List<FunctionalData> fdList) {
        List<FunctionalData> filtered = filterFunctionalData (fdList);
        List<FunctionalData> sense = filtered.stream().filter(FunctionalData::isSense).toList();
        if (sense.size()==1) {
            return sense.getFirst();
        }
        return null;
    }

    /**
     * For purposes of validation, simplify all UTR to coding.
     * @param fd
     * @return
     */
    private FunctionalData simplifyCoding (FunctionalData fd) {
        if (fd.getLocusFunction()!=LocusFunction.UTR)
            return fd;
        FunctionalData result = new FunctionalData(fd.getGene(), fd.getGeneStrand(), LocusFunction.CODING, fd.getReadStrand());
        return result;
    }

    /**
     * During expression analysis, the raw sets of functional data are filtered to
     * those that are accepted as counts by a particular priority.
     * This can result in a few outputs:
     * 1) A list with a single element is an unambigous valid assignment
     * 2) A list with no elements does not contain a gene that can be assigned.  This happens when
     * there are no functional annotations that are valid given the filters.  This can happen for example
     * when the highest priority annotation is a coding antisense annotation, and that annotation does
     * not contribute to the counts matrix.  This can also happen when the
     * priority restricts the functional annotation to a subset like (CODING,UTR) and the highest priority
     * annotation is INTRONIC, which is then filtered.
     * 3) A list with multiple elements, where multiple genes pass filters, so the gene to use is
     * ambiguous.
     * @param fdList A list of functional data to possibly reduce.
     * @return A list of functional data that is valid to count as expression.
     */
    private List<FunctionalData> filterFunctionalData (List<FunctionalData> fdList) {
        // the processes that happen during normal expression.
        // This filters to accepted annotations and simplifies compound genes.
        List<FunctionalData> filtered = gfp.getFdp().getFilteredFunctionalData(fdList);
        // this filters on strand, so antisense functional annotations are removed.
        filtered = gfp.getFdp().filterToPreferredAnnotations(filtered);
        return filtered;
    }

    private List<FunctionalData> filterFunctionalDataWithoutStrand(List<FunctionalData> fdList) {
        List<FunctionalData> filtered = gfp.getFdp().getFilteredFunctionalData(fdList);
        List<FunctionalData> result = util.filterToPreferredAnnotations(fdList, gfp.getFdp().getPriority());
        return result;
    }

    boolean containsFunctionalAnnotation(List<FunctionalData> fdList, LocusFunction func) {
        return fdList.stream().anyMatch(fd -> fd.getLocusFunction() == func);
    }


}
