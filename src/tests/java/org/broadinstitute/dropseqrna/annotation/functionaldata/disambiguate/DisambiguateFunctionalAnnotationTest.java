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
     * Add an additional functional annotation to the sense intronic read that it actually maps to another sense intronic gene.
     * We wouldn't know what gene to assign this to if the UMI resolved to being intronic, but we can unambiguously assign this to
     * an antisense coding region.
     * There's an additional flag for this being a complex region because of some strange mapping, so these
     * read sets can be filtered out of further analysis.
     */
    @Test
    public void testComplexAntiSenseCoding() {
        Map<String, List<FunctionalData>> recs = new HashMap<>();

        // the ambiguous read that maps to intronic sense and coding antisense + intronic antisense
        // but there's two intronic genes
        List<FunctionalData> ambiguousRead = FunctionalData.buildFD(
                new String [] {"A", "B", "C"},
                new String [] {"-", "+", "+"},
                new LocusFunction[] {LocusFunction.CODING, LocusFunction.INTRONIC, LocusFunction.INTRONIC},
                strandStrategy, Arrays.asList(acceptedFunctions), false);
        recs.put("complex_ambiguous", ambiguousRead);

        // the unambiguous read that maps to coding antisense only.
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
        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("C"),1);
        Assert.assertEquals(score.getAntisenseCodingCount().getCountForKey("A"),1);
        Assert.assertEquals(score.getSenseIntronicCount().getCountForKey("B"),0);
        Assert.assertTrue(score.hasAmbiguousReads());
        Assert.assertEquals(score.classify(), FunctionCategory.ANTISENSE_CODING);
        Assert.assertTrue(score.isComplex());
    }

    @Test
    /**
     * This is the mirror of the testComplexAntiSenseCoding test.
     * The resolving read is sense intronic, and even though the unambiguous read maps to
     * multiple genes, they are all sense intronic.
     */
    public void testComplexIntronicCoding() {
        Map<String, List<FunctionalData>> recs = new HashMap<>();

        // the ambiguous read that maps to intronic sense and coding antisense + intronic antisense
        // but there's two intronic genes
        List<FunctionalData> ambiguousRead = FunctionalData.buildFD(
                new String [] {"A", "B", "C"},
                new String [] {"-", "+", "+"},
                new LocusFunction[] {LocusFunction.CODING, LocusFunction.INTRONIC, LocusFunction.INTRONIC},
                strandStrategy, Arrays.asList(acceptedFunctions), false);
        recs.put("complex_ambiguous", ambiguousRead);

        // the unambiguous read that maps to the intronic sense genes B and C.
        // we don't care which gene it is.
        List<FunctionalData> unambiguousRead = FunctionalData.buildFD(
                new String [] {"B", "C"},
                new String [] {"+", "+"},
                new LocusFunction[] {LocusFunction.INTRONIC, LocusFunction.INTRONIC},
                strandStrategy, Arrays.asList(acceptedFunctions), false);
        recs.put("unambiguous", unambiguousRead);

        DisambiguateFunctionalAnnotation dfa = new DisambiguateFunctionalAnnotation();
        DisambiguationScore score= dfa.run("fake", "fake", recs);

        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("A"),1);
        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("B"),1);
        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("C"),1);
        Assert.assertEquals(score.getAntisenseCodingCount().getCountForKey("A"),0);
        Assert.assertEquals(score.getSenseIntronicCount().getCountForKey("B"),1);
        Assert.assertEquals(score.getSenseIntronicCount().getCountForKey("C"),1);
        Assert.assertTrue(score.hasAmbiguousReads());
        Assert.assertEquals(score.classify(), FunctionCategory.SENSE_INTRONIC);
        Assert.assertTrue(score.isComplex());
    }

    /**
     * This is a case that looks complex, but the complex read is actually a simple
     * coding interpretation, even without the "unambiguous" read.
     */
    @Test
    public void testComplexUnambiguous() {
        Map<String, List<FunctionalData>> recs = new HashMap<>();

        // the ambiguous read that maps to intronic sense and coding antisense + intronic antisense
        // but there's two second intronic genes, which makes the whole thing ambiguous.
        List<FunctionalData> ambiguousRead = FunctionalData.buildFD(
                new String [] {"A", "B", "C"},
                new String [] {"-", "+", "+"},
                new LocusFunction[] {LocusFunction.CODING, LocusFunction.INTRONIC, LocusFunction.CODING},
                strandStrategy, Arrays.asList(acceptedFunctions), false);
        recs.put("complex_ambiguous", ambiguousRead);

        DisambiguateFunctionalAnnotation dfa = new DisambiguateFunctionalAnnotation();
        DisambiguationScore score= dfa.run("fake", "fake", recs);

        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("A"),0);
        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("B"),0);
        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("C"),0);
        Assert.assertEquals(score.getAntisenseCodingCount().getCountForKey("A"),0);
        Assert.assertEquals(score.getSenseIntronicCount().getCountForKey("B"),0);
        Assert.assertFalse(score.hasAmbiguousReads());
        Assert.assertEquals(score.classify(), FunctionCategory.UNAMBIGUOUS);
        Assert.assertFalse(score.isComplex());

        // add an unamiguous read.  The result doesn't change.
        List<FunctionalData> unambiguousRead = FunctionalData.buildFD(
                new String [] {"C"},
                new String [] {"+"},
                new LocusFunction[] {LocusFunction.CODING},
                strandStrategy, Arrays.asList(acceptedFunctions), false);
        recs.put("unambiguous", unambiguousRead);
        score= dfa.run("fake", "fake", recs);
        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("A"),0);
        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("B"),0);
        Assert.assertEquals(score.getAmbiguousCount().getCountForKey("C"),0);
        Assert.assertEquals(score.getAntisenseCodingCount().getCountForKey("A"),0);
        Assert.assertEquals(score.getSenseIntronicCount().getCountForKey("B"),0);
        Assert.assertFalse(score.hasAmbiguousReads());
        Assert.assertEquals(score.classify(), FunctionCategory.UNAMBIGUOUS);
        Assert.assertFalse(score.isComplex());

    }


}