package org.broadinstitute.dropseqrna.annotation.functionaldata;

import org.testng.Assert;
import org.testng.annotations.Test;
import picard.annotation.LocusFunction;

public class FunctionalDataTest {

    @Test public void testType() {
        FunctionalData senseCoding = new FunctionalData("A", "+", LocusFunction.CODING, "+");
        Assert.assertEquals(FunctionalData.Type.CODING_SENSE, senseCoding.getType());

        FunctionalData senseCoding2 = new FunctionalData("A", "+", LocusFunction.UTR, "+");
        Assert.assertEquals(FunctionalData.Type.CODING_SENSE, senseCoding2.getType());

        FunctionalData antisenseCoding = new FunctionalData("A", "+", LocusFunction.CODING, "-");
        Assert.assertEquals(FunctionalData.Type.CODING_ANTISENSE, antisenseCoding.getType());

        FunctionalData antisenseCoding2 = new FunctionalData("A", "+", LocusFunction.UTR, "-");
        Assert.assertEquals(FunctionalData.Type.CODING_ANTISENSE, antisenseCoding2.getType());

        FunctionalData senseIntronic = new FunctionalData("A", "+", LocusFunction.INTRONIC, "+");
        Assert.assertEquals(FunctionalData.Type.INTRONIC_SENSE, senseIntronic.getType());

        FunctionalData antisenseIntronic = new FunctionalData("A", "+", LocusFunction.INTRONIC, "-");
        Assert.assertEquals(FunctionalData.Type.INTRONIC_ANTISENSE, antisenseIntronic.getType());

        // lots of intergenic options - strand doesn't matter, "catch all" for other classifications
        FunctionalData intergenic1 = new FunctionalData("A", "+", LocusFunction.RIBOSOMAL, "-");
        Assert.assertEquals(FunctionalData.Type.INTERGENIC, intergenic1.getType());

        FunctionalData intergenic2 = new FunctionalData("A", "+", LocusFunction.RIBOSOMAL, "+");
        Assert.assertEquals(FunctionalData.Type.INTERGENIC, intergenic2.getType());

        FunctionalData intergenic3 = new FunctionalData("A", "+", LocusFunction.INTERGENIC, "-");
        Assert.assertEquals(FunctionalData.Type.INTERGENIC, intergenic3.getType());

        FunctionalData intergenic4 = new FunctionalData("A", "+", LocusFunction.INTERGENIC, "+");
        Assert.assertEquals(FunctionalData.Type.INTERGENIC, intergenic4.getType());

    }
}