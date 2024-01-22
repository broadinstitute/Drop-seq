package org.broadinstitute.dropseqrna.annotation.functionaldata.disambiguate;

import org.broadinstitute.dropseqrna.annotation.functionaldata.*;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import picard.annotation.LocusFunction;

import java.util.*;
import java.util.stream.Collectors;

public class DisambiguateFunctionalAnnotation {

    private static final List<LocusFunction> DEFAULT_LOCUS_FUNCTION_LIST = List.of(LocusFunction.CODING, LocusFunction.UTR, LocusFunction.INTRONIC);
    private static final StrandStrategy DEFAULT_STRAND_STRATEGY = StrandStrategy.SENSE;

    private final DataProcessorUtils util;

    private final StarSoloPriorityScore priorityScore = new StarSoloPriorityScore();

    /**
     * The default test, which compares the DropSeq and StarSolo strategies.
     * This considers Coding+Intronic reads on the sense strand.
     */
    public DisambiguateFunctionalAnnotation () {
        util = new DataProcessorUtils(DEFAULT_STRAND_STRATEGY, DEFAULT_LOCUS_FUNCTION_LIST);
    }
    /**
     * @param recs A list of SAMRecord names for a single cell and UMI mapped to their functional annotations
     * @return Return the structure containing the counts of ambiguous and unambiguous classified reads, one for each input read.
     */
    public DisambiguationScore run (String cellBarcode, String molecularBarcode, Map<String, List<FunctionalData>> recs) {
        // TODO: fast fail if there is only a single record for this cell/umi.
        DisambiguationScore score = new DisambiguationScore(cellBarcode, molecularBarcode);

        // StarSolo should emit antisense functional annotations.
        // StarSolo alone should provide the right set of interpretations.
        Map<String, List<FunctionalData>> fdMap = simplify(recs);

        // walk through the reads and look for cases where reads are ambiguous or support an unambiguous call.
        for (String readName: fdMap.keySet()) {
            List<FunctionalData> fd = fdMap.get(readName);
            // score the list of functional data for this read.
            score.add(fd);
        }
        return score;
    }

    public Map<String, List<FunctionalData>>  simplify (Map<String, List<FunctionalData>> recs) {

        Map<String, List<FunctionalData>> result = new HashMap<>();
        for (String readName: recs.keySet()) {
            List<FunctionalData> fdList = recs.get(readName);
            // filter on accepted functions
            fdList = util.filterOnLocus(fdList);
            // simplify those functions but in a strand-specific way
            fdList = filterToPreferredAnnotationsStrandSpecific(fdList);
            result.put(readName, fdList);
        }
        return result;
    }

    private static <T> boolean haveEqualSets(List<T> list1, List<T> list2) {
        Set<T> set1 = new HashSet<>(list1);
        Set<T> set2 = new HashSet<>(list2);

        return set1.equals(set2);
    }


    /**
     * For each strand, simplify functional annotations to the preferred one.
     * This simplifies CODING+INTRONIC -> CODING.
     * But if coding and intronic occur on opposite strands, they are not simplified.
     * @return The list of functional annotations, simplified to the highest priority annotation(s) per strand.
     */
//    private List<FunctionalData> filterToPreferredAnnotationsStrandSpecificOld (List<FunctionalData> fdList) {
//        List<FunctionalData> result = new ArrayList<>();
//
//        Map<String, List<FunctionalData>> strandMap = fdList.stream()
//                .collect(Collectors.groupingBy(FunctionalData::getGeneStrand));
//
//        for (String strand: strandMap.keySet()) {
//            List<FunctionalData> fd = strandMap.get(strand);
//            fd=util.filterToPreferredAnnotations(fd, priorityScore);
//            result.addAll(fd);
//        }
//        return result;
//    }

    /**
     * For each strand, simplify functional annotations to the preferred one.
     * This simplifies CODING+INTRONIC -> CODING.
     * But if coding and intronic occur on opposite strands, they are not simplified.
     *
     * Thanks chatgpt!
     *
     * @return The list of functional annotations, simplified to the highest priority annotation(s) per strand.
     */    private List<FunctionalData> filterToPreferredAnnotationsStrandSpecific(List<FunctionalData> fdList) {
        return fdList.stream()
                .collect(Collectors.groupingBy(FunctionalData::getGeneStrand))
                .values().stream()
                .flatMap(strandGroup -> util.filterToPreferredAnnotations(strandGroup, priorityScore).stream())
                .collect(Collectors.toList());
    }


    boolean isAntisenseCoding (FunctionalData d) {
        return (d.getLocusFunction()==LocusFunction.CODING |d.getLocusFunction()==LocusFunction.UTR) & !d.getGeneStrand().equals(d.getReadStrand());
    }

    boolean isSenseIntronic (FunctionalData d) {
        return d.getLocusFunction()==LocusFunction.INTRONIC & d.getGeneStrand()==d.getReadStrand();
    }



}
