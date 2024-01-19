package org.broadinstitute.dropseqrna.annotation.functionaldata.disambiguate;

import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalData;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.annotation.LocusFunction;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class DisambiguateFunctionalAnnotationTest {

    private final LocusFunction [] acceptedFunctions = {LocusFunction.CODING, LocusFunction.UTR, LocusFunction.INTRONIC};
    private final StrandStrategy strandStrategy = StrandStrategy.SENSE;

    @Test
    /**
     * No ambiguous reads.  One coding read, one intronic read.
     */
    public void testUnambiguous() {
        Map<String, List<FunctionalData>> recs = new HashMap<>();

        // 2 reads where the gene and read strand match.
        List<FunctionalData> unambiguousRead = FunctionalData.buildFD(
                new String [] {"A"},
                new String [] {"-"},
                new LocusFunction[] {LocusFunction.CODING},
                strandStrategy, Arrays.asList(acceptedFunctions), true);
        recs.put("unambiguous", unambiguousRead);

        List<FunctionalData> unambiguousRead2 = FunctionalData.buildFD(
                new String [] {"B"},
                new String [] {"-"},
                new LocusFunction[] {LocusFunction.INTRONIC},
                strandStrategy, Arrays.asList(acceptedFunctions), true);
        recs.put("unambiguous2", unambiguousRead2);

        DisambiguateFunctionalAnnotation dfa = new DisambiguateFunctionalAnnotation();
        DisambiguationScore score= dfa.run("fake", "fake", recs);

        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("A"),0);
        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("B"),0);
        Assert.assertEquals(score.getSenseIntronicCount().getCountForKey("B"),1);

        Assert.assertFalse(score.hasAmbiguousReads());
        Assert.assertEquals(score.classify(), FunctionCategory.UNAMBIGUOUS);
        Assert.assertFalse(score.isComplex());
    }

    @Test
    /**
     * Simple case of a single ambiguous read and a single antisense coding read that is successfully disambiguated.
     */
    public void testBasicDisambiguation() {
        Map<String, List<FunctionalData>> recs = new HashMap<>();

        // the ambiguous read that maps to intronic sense and coding antisense
        List<FunctionalData> ambiguousRead = FunctionalData.buildFD(
                new String [] {"A", "B"},
                new String [] {"-", "+"},
                new LocusFunction[] {LocusFunction.CODING, LocusFunction.INTRONIC},
                strandStrategy, Arrays.asList(acceptedFunctions), false);

        // the unambiguous read that maps to coding antisense only.  This would resolve the UMI.
        List<FunctionalData> unambiguousRead = FunctionalData.buildFD(
                new String [] {"A"},
                new String [] {"-"},
                new LocusFunction[] {LocusFunction.CODING},
                strandStrategy, Arrays.asList(acceptedFunctions), false);

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
        Assert.assertFalse(score.isComplex());
    }

    @Test
    /**
     * Simple case of a single ambiguous read and a single antisense coding read that is successfully disambiguated.
     * Add in unambiguous intronic and antisense coding reads on other genes that do not impact the ambiguous gene
     */
    public void testNoisyBasicDisambiguation() {
        Map<String, List<FunctionalData>> recs = new HashMap<>();

        // the ambiguous read that maps to intronic sense and coding antisense
        List<FunctionalData> ambiguousRead = FunctionalData.buildFD(
                new String [] {"A", "B"},
                new String [] {"-", "+"},
                new LocusFunction[] {LocusFunction.CODING, LocusFunction.INTRONIC},
                strandStrategy, Arrays.asList(acceptedFunctions), false);
        recs.put("ambiguous", ambiguousRead);
        // the unambiguous read that maps to coding antisense only.  This would resolve the UMI.
        List<FunctionalData> unambiguousRead = FunctionalData.buildFD(
                new String [] {"A"},
                new String [] {"-"},
                new LocusFunction[] {LocusFunction.CODING},
                strandStrategy, Arrays.asList(acceptedFunctions), false);
        recs.put("unambiguous", unambiguousRead);

        // add some junky noise in - an intronic read on a different gene (C)
        // and an antisense coding read on a different gene (D)
        // these should not impact the ambiguous read of if it's disambiguated successfully

        List<FunctionalData> noiseSenseIntronic = FunctionalData.buildFD(
                new String [] {"C"},
                new String [] {"+"},
                new LocusFunction[] {LocusFunction.INTRONIC},
                strandStrategy, Arrays.asList(acceptedFunctions), false);
        recs.put("noiseSenseIntronic", noiseSenseIntronic);

        List<FunctionalData> noiseAntisenseCoding = FunctionalData.buildFD(
                new String [] {"D"},
                new String [] {"-"},
                new LocusFunction[] {LocusFunction.CODING},
                strandStrategy, Arrays.asList(acceptedFunctions), false);
        recs.put("noiseAntisenseCoding", noiseAntisenseCoding);


        DisambiguateFunctionalAnnotation dfa = new DisambiguateFunctionalAnnotation();
        DisambiguationScore score= dfa.run("fake", "fake", recs);

        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("A"),1);
        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("B"),1);
        Assert.assertEquals(score.getAntisenseCodingCount().getCountForKey("A"),1);
        Assert.assertEquals(score.getSenseIntronicCount().getCountForKey("B"),0);
        Assert.assertEquals(score.getAntisenseCodingCount().getCountForKey("D"),1);
        Assert.assertEquals(score.getSenseIntronicCount().getCountForKey("C"),1);
        Assert.assertTrue(score.hasAmbiguousReads());
        Assert.assertEquals(score.classify(), FunctionCategory.ANTISENSE_CODING);
        Assert.assertFalse(score.isComplex());
    }


    /**
     * Set up an ambiguous gene pair.
     * Add additional observations to the antisense coding read (antisense intronic) that needs to be "simplified"
     * from coding+introic to coding first before testing.
     */
    @Test
    public void testSimplificationDisambiguation() {
        Map<String, List<FunctionalData>> recs = new HashMap<>();

        // the ambiguous read that maps to intronic sense and coding antisense + intronic antisense
        List<FunctionalData> ambiguousRead = FunctionalData.buildFD(
                new String [] {"A", "B", "A"},
                new String [] {"-", "+", "-"},
                new LocusFunction[] {LocusFunction.CODING, LocusFunction.INTRONIC, LocusFunction.INTRONIC},
                strandStrategy, Arrays.asList(acceptedFunctions), false);
        recs.put("ambiguous", ambiguousRead);
        // the unambiguous read that maps to coding antisense only.  This would resolve the UMI.
        List<FunctionalData> unambiguousRead = FunctionalData.buildFD(
                new String [] {"A"},
                new String [] {"-"},
                new LocusFunction[] {LocusFunction.CODING},
                strandStrategy, Arrays.asList(acceptedFunctions), false);
        recs.put("unambiguous", unambiguousRead);

        DisambiguateFunctionalAnnotation dfa = new DisambiguateFunctionalAnnotation();
        DisambiguationScore score= dfa.run("fake", "fake", recs);

        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("A"),1);
        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("B"),1);
        Assert.assertEquals(score.getAntisenseCodingCount().getCountForKey("A"),1);
        Assert.assertEquals(score.getSenseIntronicCount().getCountForKey("B"),0);
        Assert.assertTrue(score.hasAmbiguousReads());
        Assert.assertEquals(score.classify(), FunctionCategory.ANTISENSE_CODING);
        Assert.assertFalse(score.isComplex());
    }

    /**
     * Set up an ambiguous gene pair.
     * Add an additional functional annotation to the sense intronic read that it actually maps to another sense intronic gene,
     * and thus is uninformative because it's unclear which of the two intronic reads is being used.  This should NOT resolve the ambiguity.
     * There's an additional flag for this being a complex region because of some strange mapping, so these
     * read sets can be filtered out of further analysis.
     */
    @Test
    public void testComplexAmbiguous1() {
        Map<String, List<FunctionalData>> recs = new HashMap<>();

        // the ambiguous read that maps to intronic sense and coding antisense + intronic antisense
        // but there's two second intronic genes, which makes the whole thing ambiguous.
        List<FunctionalData> ambiguousRead = FunctionalData.buildFD(
                new String [] {"A", "B", "C"},
                new String [] {"-", "+", "+"},
                new LocusFunction[] {LocusFunction.CODING, LocusFunction.INTRONIC, LocusFunction.INTRONIC},
                strandStrategy, Arrays.asList(acceptedFunctions), false);
        recs.put("complex_ambiguous", ambiguousRead);

        // the unambiguous read that maps to coding antisense only.  This would resolve the UMI if
        // the ambiguous read wasn't complex.
        List<FunctionalData> unambiguousRead = FunctionalData.buildFD(
                new String [] {"A"},
                new String [] {"-"},
                new LocusFunction[] {LocusFunction.CODING},
                strandStrategy, Arrays.asList(acceptedFunctions), false);
        recs.put("unambiguous", unambiguousRead);

        DisambiguateFunctionalAnnotation dfa = new DisambiguateFunctionalAnnotation();
        DisambiguationScore score= dfa.run("fake", "fake", recs);

        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("A"),0);
        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("B"),0);
        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("C"),0);
        Assert.assertEquals(score.getAntisenseCodingCount().getCountForKey("A"),1);
        Assert.assertEquals(score.getSenseIntronicCount().getCountForKey("B"),0);
        Assert.assertFalse(score.hasAmbiguousReads()); // can't disambiguate B and C.
        Assert.assertEquals(score.classify(), FunctionCategory.UNAMBIGUOUS);
        Assert.assertTrue(score.isComplex());
    }

    /**
     * Set up an ambiguous gene pair.
     * Add an additional functional annotation to the sense intronic read that maps to a sense coding gene.
     * This means the comparison should be between a coding gene on the correct strand and a coding gene on the incorrect strand.
     * Thus, this should not be ambiguous in the first place, the coding gene should be used by DGE.
    */
    @Test
    public void testComplexUnambiguous1() {
        Map<String, List<FunctionalData>> recs = new HashMap<>();

        // the ambiguous read that maps to intronic sense and coding antisense + intronic antisense
        // but there's two second intronic genes, which makes the whole thing ambiguous.
        List<FunctionalData> ambiguousRead = FunctionalData.buildFD(
                new String [] {"A", "B", "B"},
                new String [] {"-", "+", "+"},
                new LocusFunction[] {LocusFunction.CODING, LocusFunction.INTRONIC, LocusFunction.CODING},
                strandStrategy, Arrays.asList(acceptedFunctions), false);
        recs.put("complex_ambiguous", ambiguousRead);

        // the unambiguous read that maps to coding antisense only.  This would resolve the UMI if
        // the ambiguous read wasn't complex.
        List<FunctionalData> unambiguousRead = FunctionalData.buildFD(
                new String [] {"A"},
                new String [] {"-"},
                new LocusFunction[] {LocusFunction.CODING},
                strandStrategy, Arrays.asList(acceptedFunctions), false);
        recs.put("unambiguous", unambiguousRead);

        DisambiguateFunctionalAnnotation dfa = new DisambiguateFunctionalAnnotation();
        DisambiguationScore score= dfa.run("fake", "fake", recs);

        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("A"),0);
        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("B"),0);
        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("C"),0);
        Assert.assertEquals(score.getAntisenseCodingCount().getCountForKey("A"),1);
        Assert.assertEquals(score.getSenseIntronicCount().getCountForKey("B"),0);
        Assert.assertFalse(score.hasAmbiguousReads()); // can't disambiguate B and C.
        Assert.assertEquals(score.classify(), FunctionCategory.UNAMBIGUOUS);
        Assert.assertFalse(score.isComplex());
    }


}