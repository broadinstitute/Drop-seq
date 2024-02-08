package org.broadinstitute.dropseqrna.annotation.functionaldata.disambiguate;

import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalData;
import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalDataProcessorStrategy;
import org.broadinstitute.dropseqrna.barnyard.GeneFunctionCommandLineBase;
import org.broadinstitute.dropseqrna.utils.readiterators.GeneFunctionProcessor;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.annotation.LocusFunction;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class ClassifyDropSeqFunctionalDataTest {

    private static String NO_GENE = ClassifyDropSeqFunctionalData.DEFAULT_MISSING_GENE_LABEL;
    private ClassifyDropSeqFunctionalData getDefaultValidator () {
        List<LocusFunction> LOCUS_FUNCTION_LIST = Collections.unmodifiableList(new ArrayList<>(Arrays.asList(LocusFunction.CODING, LocusFunction.UTR, LocusFunction.INTRONIC)));
        StrandStrategy STRAND_STRATEGY=GeneFunctionCommandLineBase.DEFAULT_STRAND_STRATEGY;

        GeneFunctionProcessor gfp = new GeneFunctionProcessor(GeneFunctionCommandLineBase.DEFAULT_GENE_NAME_TAG, GeneFunctionCommandLineBase.DEFAULT_GENE_STRAND_TAG,
                GeneFunctionCommandLineBase.DEFAULT_GENE_FUNCTION_TAG,false,
                STRAND_STRATEGY, LOCUS_FUNCTION_LIST, FunctionalDataProcessorStrategy.STARSOLO);

        return (new ClassifyDropSeqFunctionalData(gfp));
    }



    @Test
    public void testIsDropSeqAmbiguous() {
        ClassifyDropSeqFunctionalData va = getDefaultValidator();

        // two coding on the same strand - definitely ambiguous.
        List<FunctionalData> testAmbiguousSenseCoding= constructAmbiguousSenseCoding();
        Assert.assertTrue(va.isDropSeqAmbiguous(testAmbiguousSenseCoding));

        // test the full converter
        GeneWithFunction expected = new GeneWithFunction(ClassifyDropSeqFunctionalData.DEFAULT_MISSING_GENE_LABEL, FunctionalData.Type.AMBIGUOUS);
        GeneWithFunction gf =  va.convert(testAmbiguousSenseCoding);
        Assert.assertEquals(gf, expected);

        // two intronic
        List<FunctionalData> testAmbiguousSenseIntronic= constructAmbiguousSenseIntronic();
        Assert.assertTrue(va.isDropSeqAmbiguous(testAmbiguousSenseIntronic));
        // test the full converter
        GeneWithFunction gf2 =  va.convert(testAmbiguousSenseIntronic);
        Assert.assertEquals(gf2, expected);


        // sense coding A + sense intronic B is not ambiguous.
        FunctionalData senseCodingA = new FunctionalData("A", "+", LocusFunction.CODING, "+");
        FunctionalData senseIntronicB = new FunctionalData("B", "+", LocusFunction.INTRONIC, "+");
        List<FunctionalData> testAmbiguous2 = Arrays.asList(senseCodingA, senseIntronicB);
        Assert.assertFalse(va.isDropSeqAmbiguous(testAmbiguous2));
        // this should classify as coding A.
        expected = new GeneWithFunction("A", FunctionalData.Type.CODING_SENSE);
        GeneWithFunction codingA =  va.convert(testAmbiguous2);
        Assert.assertEquals(codingA, expected);


        // sense coding A + sense intronic A is not ambiguous.
        FunctionalData senseIntronicA = new FunctionalData("A", "+", LocusFunction.INTRONIC, "+");
        List<FunctionalData> testAmbiguous3 = Arrays.asList(senseCodingA, senseIntronicA);
        Assert.assertFalse(va.isDropSeqAmbiguous(testAmbiguous3));
        // same expectation - coding A.
        GeneWithFunction codingA2 =  va.convert(testAmbiguous3);
        Assert.assertEquals(codingA2, expected);


        // antisense encoding, should not be ambiguous.
        FunctionalData antisenseCodingA = new FunctionalData("A", "+", LocusFunction.CODING, "-");
        List<FunctionalData> antisenseCoding = Arrays.asList(antisenseCodingA);
        Assert.assertFalse(va.isDropSeqAmbiguous(antisenseCoding));
        // this should classify as antisense A.
        expected = new GeneWithFunction("A", FunctionalData.Type.CODING_ANTISENSE);
        GeneWithFunction antiSenseA =  va.convert(antisenseCoding);
        Assert.assertEquals(antiSenseA, expected);

    }

    @Test
    public void testIsDropSeqAntisenseCoding() {
        ClassifyDropSeqFunctionalData va = getDefaultValidator();
        // two on the same strand and coding
        List<FunctionalData> antiSenseCoding= constructAntisenseCoding();
        Assert.assertNotNull(va.getDropSeqAntisenseCoding(antiSenseCoding));

        // test the full converter
        GeneWithFunction expected = new GeneWithFunction("A", FunctionalData.Type.CODING_ANTISENSE);
        GeneWithFunction gf =  va.convert(antiSenseCoding);
        Assert.assertEquals(gf, expected);

        // two on the same strand and intronic
        List<FunctionalData> antiSenseIntronic= constructAntisenseIntronic();
        Assert.assertNotNull(va.getDropSeqAntisenseIntronic(antiSenseIntronic));

        // test the full converter
        expected = new GeneWithFunction("A", FunctionalData.Type.INTRONIC_ANTISENSE);
        gf =  va.convert(antiSenseIntronic);
        Assert.assertEquals(gf, expected);

        // differentiate from ambiguous
        List<FunctionalData> ambiguousCoding= constructAmbiguousSenseCoding();
        Assert.assertNull(va.getDropSeqAntisenseCoding(ambiguousCoding));
        Assert.assertNull(va.getDropSeqAntisenseIntronic(ambiguousCoding));

        // test the full converter
        expected = new GeneWithFunction(ClassifyDropSeqFunctionalData.DEFAULT_MISSING_GENE_LABEL, FunctionalData.Type.AMBIGUOUS);
        gf =  va.convert(ambiguousCoding);
        Assert.assertEquals(gf, expected);

        // differentiate from a sense and antisense encoding.
        // this has a valid coding record you should get instead.
        FunctionalData A = new FunctionalData("A", "+", LocusFunction.CODING, "+");
        FunctionalData B = new FunctionalData("B", "+", LocusFunction.CODING, "-");
        List<FunctionalData> testCoding = Arrays.asList(A, B);
        Assert.assertNull(va.getDropSeqAntisenseCoding(testCoding));
        // test the full converter
        expected = new GeneWithFunction("A", FunctionalData.Type.CODING_SENSE);
        gf =  va.convert(testCoding);
        Assert.assertEquals(gf, expected);


    }

    @Test
    public void testIsDropSeqAntisenseIntronic() {
        ClassifyDropSeqFunctionalData va = getDefaultValidator();
        // two on the same strand and coding
        List<FunctionalData> antiSenseCoding= constructAntisenseIntronic();
        Assert.assertNotNull(va.getDropSeqAntisenseIntronic(antiSenseCoding));

        // test the full converter
        GeneWithFunction expected = new GeneWithFunction("A", FunctionalData.Type.INTRONIC_ANTISENSE);
        GeneWithFunction gf =  va.convert(antiSenseCoding);
        Assert.assertEquals(gf, expected);
    }

    @Test
    public void testCodingSense () {
        ClassifyDropSeqFunctionalData va = getDefaultValidator();

        FunctionalData senseCodingA = new FunctionalData("A", "+", LocusFunction.CODING, "+");
        FunctionalData senseIntronicB = new FunctionalData("B", "+", LocusFunction.INTRONIC, "+");
        List<FunctionalData> senseCoding = Arrays.asList(senseCodingA, senseIntronicB);
        Assert.assertNotNull(va.getDropSeqSense(senseCoding));

        GeneWithFunction expected = new GeneWithFunction("A", FunctionalData.Type.CODING_SENSE);
        GeneWithFunction gf =  va.convert(senseCoding);
        Assert.assertEquals(gf, expected);

        // the antisense intronic is discarded.
        FunctionalData antisenseIntronicB = new FunctionalData("B", "+", LocusFunction.INTRONIC, "-");
        List<FunctionalData> senseCoding2 = Arrays.asList(senseCodingA, antisenseIntronicB);
        Assert.assertNotNull(va.getDropSeqSense(senseCoding));

        expected = new GeneWithFunction("A", FunctionalData.Type.CODING_SENSE);
        gf =  va.convert(senseCoding);
        Assert.assertEquals(gf, expected);

        // mix in an intergenic annotation for the same gene.  This should resolve to intergenic, NOT coding if it's the same gene.
        FunctionalData intergenicA = new FunctionalData("A", "+", LocusFunction.INTERGENIC, "+");
        List<FunctionalData> intergeic = Arrays.asList(senseCodingA, intergenicA);
        Assert.assertNull(va.getDropSeqSense(intergeic));

        expected = new GeneWithFunction(ClassifyDropSeqFunctionalData.DEFAULT_MISSING_GENE_LABEL, FunctionalData.Type.INTERGENIC);
        gf =  va.convert(intergeic);
        Assert.assertEquals(gf, expected);

    }



    // a read that maps to antisense intronic twice.
    // Star for some reason tags this as antisense intronic and gives it a gene name of '-' so it's not counted.
    // LH00118:69:22CNFNLT3:1:1255:24097:2154	0	chr2	59347496	255	90M	*	0	0  CB:Z:AAACCCACAATGTCTG	UB:Z:AAATACGGTCTA	XF:Z:INTERGENIC	sF:B:i,6,0	NH:i:1
    // GN:Z:-	CR:Z:AAACCCACAATGTCTG	UR:Z:AAATACGGTCTA GX:Z:-	CY:Z:FFFFFFFFFFFFFFFF	UY:Z:FFFFFFFFFFFF	gf:Z:INTRONIC,INTRONIC	gn:Z:AC007100.1,AC007179.2	gs:Z:-,-
    @Test
    public void testAmbiguousAntisenseIntronic () {
        ClassifyDropSeqFunctionalData va = getDefaultValidator();
        FunctionalData A = new FunctionalData("A", "+", LocusFunction.INTRONIC, "-");
        FunctionalData B = new FunctionalData("B", "+", LocusFunction.INTRONIC, "-");
        List<FunctionalData> ambiguous = Arrays.asList(A, B);
        Assert.assertNull(va.getDropSeqAntisenseIntronic(ambiguous));
        Assert.assertTrue(va.isDropSeqAmbiguous(ambiguous));

        GeneWithFunction expected = new GeneWithFunction(ClassifyDropSeqFunctionalData.DEFAULT_MISSING_GENE_LABEL, FunctionalData.Type.AMBIGUOUS);
        GeneWithFunction gf =  va.convert(ambiguous);
        Assert.assertEquals(gf, expected);

    }

    @Test
    public void testCodingPlusUTRAntisense () {
        ClassifyDropSeqFunctionalData va = getDefaultValidator();
        FunctionalData a1 = new FunctionalData("A", "+", LocusFunction.CODING, "-");
        FunctionalData a2 = new FunctionalData("A", "+", LocusFunction.UTR, "-");

        List<FunctionalData> antisense = Arrays.asList(a1, a2);
        Assert.assertNotNull(va.getDropSeqAntisenseCoding(antisense));
        Assert.assertFalse(va.isDropSeqAmbiguous(antisense));

        GeneWithFunction expected = new GeneWithFunction("A", FunctionalData.Type.CODING_ANTISENSE);
        GeneWithFunction gf =  va.convert(antisense);
        Assert.assertEquals(gf, expected);

    }

    @Test
    /**
     * Classified as ambiguous / intronic.
     * The coding overlaps intergenic so should be removed, leaving two genes that are intronic sense.
     * 0 = {FunctionalData@3570} "AL122008.1 + INTRONIC + INTRONIC_SENSE"
     * 1 = {FunctionalData@3571} "AL122008.4 - CODING + CODING_ANTISENSE"
     * 2 = {FunctionalData@3572} "AL122008.4 - INTERGENIC + INTERGENIC"
     * 3 = {FunctionalData@3573} "SLC35F3 + INTRONIC + INTRONIC_SENSE"
     */
    public void complexAntisense () {
        ClassifyDropSeqFunctionalData va = getDefaultValidator();
        FunctionalData a = new FunctionalData("A", "+", LocusFunction.INTRONIC, "+");
        FunctionalData b1= new FunctionalData("B", "+", LocusFunction.CODING, "-");
        FunctionalData b2= new FunctionalData("B", "+", LocusFunction.INTERGENIC, "-");
        FunctionalData c= new FunctionalData("C", "+", LocusFunction.INTRONIC, "+");

        List<FunctionalData> test = Arrays.asList(a, b1,b2,c);
        Assert.assertNull(va.getDropSeqAntisenseCoding(test));
        Assert.assertTrue(va.isDropSeqAmbiguous(test));
        Assert.assertNull(va.getDropSeqSense(test));

        GeneWithFunction expected = new GeneWithFunction(NO_GENE, FunctionalData.Type.AMBIGUOUS);
        GeneWithFunction gf =  va.convert(test);
        Assert.assertEquals(gf, expected);

    }



    private List<FunctionalData> constructAntisenseCoding () {
        // ANTISENSE CODING + SENSE INTRONIC -> antisense when using STARSolo rules.
        FunctionalData senseCodingA = new FunctionalData("A", "+", LocusFunction.CODING, "-");
        FunctionalData senseCodingB = new FunctionalData("B", "+", LocusFunction.INTRONIC, "+");
        List<FunctionalData> antiSenseCoding = Arrays.asList(senseCodingA, senseCodingB);
        return antiSenseCoding;
    }

    private List<FunctionalData> constructAntisenseIntronic () {
        FunctionalData A = new FunctionalData("A", "+", LocusFunction.INTRONIC, "-");
        List<FunctionalData> antiSenseIntronic = List.of(A);
        return antiSenseIntronic;
    }


    // two sense coding on the same strand
    private List<FunctionalData> constructAmbiguousSenseCoding () {
        // two coding on the same strand - definitely ambiguous.
        FunctionalData senseCodingA = new FunctionalData("A", "+", LocusFunction.CODING, "+");
        FunctionalData senseCodingB = new FunctionalData("B", "+", LocusFunction.CODING, "+");
        List<FunctionalData> testAmbiguous = Arrays.asList(senseCodingA, senseCodingB);
        return testAmbiguous;
    }

    // two sense intronic on the same strand
    private List<FunctionalData> constructAmbiguousSenseIntronic () {
        // two coding on the same strand - definitely ambiguous.
        FunctionalData senseCodingA = new FunctionalData("A", "+", LocusFunction.INTRONIC, "+");
        FunctionalData senseCodingB = new FunctionalData("B", "+", LocusFunction.INTRONIC, "+");
        List<FunctionalData> testAmbiguous = Arrays.asList(senseCodingA, senseCodingB);
        return testAmbiguous;
    }
}