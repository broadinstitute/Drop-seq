package org.broadinstitute.dropseqrna.annotation.functionaldata;

import org.apache.commons.collections4.map.HashedMap;
import picard.annotation.LocusFunction;

import java.util.*;

public class StarSoloPriorityScore implements PriorityScoreI {

    public int getScore(FunctionalData fd) {
        int result = switch (fd.getLocusFunction()) {
            case LocusFunction.CODING -> 1;
            case LocusFunction.UTR -> 1;
            case LocusFunction.INTRONIC -> 2;
            case LocusFunction.INTERGENIC -> 3;
            default -> Integer.MAX_VALUE;
        };

        return result;
    }
}


