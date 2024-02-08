package org.broadinstitute.dropseqrna.annotation.functionaldata;

import picard.annotation.LocusFunction;

public class DropSeqPriorityScore implements PriorityScoreI {
    public int getScore(FunctionalData fd) {
        int result = switch (fd.getLocusFunction()) {
            case null -> Integer.MAX_VALUE;
            case LocusFunction.CODING -> 1;
            case LocusFunction.UTR -> 1;
            case LocusFunction.INTRONIC -> 2;
            case LocusFunction.INTERGENIC -> 3;
            default -> Integer.MAX_VALUE;
        };

        return result;
    }
}


