package org.broadinstitute.dropseqrna.annotation.functionaldata.disambiguate;

import org.broadinstitute.dropseqrna.annotation.functionaldata.*;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import picard.annotation.LocusFunction;

import java.util.*;

public class DisambiguateFunctionalAnnotation {

    private static final List<LocusFunction> DEFAULT_LOCUS_FUNCTION_LIST = List.of(LocusFunction.CODING, LocusFunction.UTR, LocusFunction.INTRONIC);
    private static final StrandStrategy DEFAULT_STRAND_STRATEGY = StrandStrategy.SENSE;

    private final DataProcessorUtils util;

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
            // simplify those functions
            fdList = util.simplifyFD(fdList);
            result.put(readName, fdList);
        }
        return result;
    }

    boolean isAntisenseCoding (FunctionalData d) {
        return (d.getLocusFunction()==LocusFunction.CODING |d.getLocusFunction()==LocusFunction.UTR) & !d.getGeneStrand().equals(d.getReadStrand());
    }

    boolean isSenseIntronic (FunctionalData d) {
        return d.getLocusFunction()==LocusFunction.INTRONIC & d.getGeneStrand()==d.getReadStrand();
    }



}
