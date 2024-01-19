package org.broadinstitute.dropseqrna.annotation.functionaldata.disambiguate;

import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalData;
import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalDataProcessorI;
import org.broadinstitute.dropseqrna.annotation.functionaldata.StarSoloFunctionalDataProcessor;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.annotation.LocusFunction;

import java.util.*;

import static org.testng.Assert.*;

public class DisambiguateFunctionalAnnotationTest {

    private LocusFunction [] acceptedFunctions = {LocusFunction.CODING, LocusFunction.UTR, LocusFunction.INTRONIC};
    private StrandStrategy strandStrategy = StrandStrategy.SENSE;

    FunctionalDataProcessorI fdp = new StarSoloFunctionalDataProcessor(StrandStrategy.SENSE, acceptedFunctions);

    @Test
    public void testBasicDisambiguation() {
        Map<String, List<FunctionalData>> recs = new HashMap<>();

        // the ambiguous read that maps to intronic sense and coding antisense
        List<FunctionalData> ambiguousRead = FunctionalData.buildFD(
                new String [] {"A", "B"},
                new String [] {"-", "+"},
                new LocusFunction[] {LocusFunction.CODING, LocusFunction.INTRONIC},
                StrandStrategy.SENSE, Arrays.asList(acceptedFunctions), false);

        // the unambiguous read that maps to coding antisense only.  This would resolve the UMI.
        List<FunctionalData> unambiguousRead = FunctionalData.buildFD(
                new String [] {"A"},
                new String [] {"-"},
                new LocusFunction[] {LocusFunction.CODING},
                StrandStrategy.SENSE, Arrays.asList(acceptedFunctions), false);

        recs.put("ambiguous", ambiguousRead);
        recs.put("unambiguous", unambiguousRead);

        DisambiguateFunctionalAnnotation dfa = new DisambiguateFunctionalAnnotation();
        DisambiguationScore score= dfa.run("fake", "fake", recs);

        Assert.assertEquals(score.getTotalCount(),3);
        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("A"),1);
        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("B"),1);
        Assert.assertEquals(score.getAntisenseCodingCount().getCountForKey("A"),1);
        Assert.assertEquals(score.getSenseIntronicCount().getCountForKey("B"),0);

        Assert.assertTrue(score.hasAmbiguousReads());
        Assert.assertEquals(score.classify(), FunctionCategory.ANTISENSE_CODING);
    }

}