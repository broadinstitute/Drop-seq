/*
 * MIT License
 *
 * Copyright 2017 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.multisample;

import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.GenotypeType;

import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.multisample.VariantData;
import org.testng.Assert;
import org.testng.annotations.Test;

public class VariantDataTest {

	@Test
	public void testImpossibleOne () {
		char [] bases = {'A', 'A', 'T', 'T', 'T'};
		int [] qualities = {10,10,10,10,10};
		VariantData vd = new VariantData(new Interval("chr1", 1, 1), 'A', 'T', GenotypeType.HOM_REF, GenotypeType.HET, bases, qualities);

		int impossibleOne=vd.getNumImpossibleAllelesGenotypeOne();
		int impossibleTwo=vd.getNumImpossibleAllelesGenotypeTwo();
		Assert.assertEquals(impossibleOne, 3);
		Assert.assertEquals(impossibleTwo, 0);

	}

	@Test(expectedExceptions=IllegalArgumentException.class)
	// test a no-call SNP state with no missing data.
	public void testMissingDataWithException() {
		char [] bases = {'A', 'A', 'A', 'A', 'T'};
		int [] qualities = {10,10,10,10,10};
		VariantData vd = new VariantData(new Interval("chr1", 1, 1), 'A', 'T', GenotypeType.HOM_REF, GenotypeType.NO_CALL, bases, qualities, null, null, null, null);
		@SuppressWarnings("unused")
		double one = vd.getLogLikelihood(1);
	}

	@Test
	public void testMissingData1() {
		char [] bases = {'A', 'A', 'A', 'A', 'T'};
		int [] qualities = {10,10,10,10,10};
		double missingDataLikelihood =0.6;
		VariantData vd = new VariantData(new Interval("chr1", 1, 1), 'A', 'T', GenotypeType.HOM_REF, GenotypeType.NO_CALL, bases, qualities, missingDataLikelihood, null, null, null);
		// this should be sum(log10(c(0.9, 0.9,0.9,0.9,0.1)))=-1.18303
		double one = vd.getLogLikelihood(1);
		Assert.assertEquals(one, -1.18303, 0.01);
		// this should be sum(log10(c(0.6,0.6,0.6,0.6,0.6)))=-1.109244
		double zero = vd.getLogLikelihood(0);
		Assert.assertEquals(zero, -1.109244, 0.01);
	}

	//TODO: DO this!
	public void testMissingData() {
		//1:43167815-43167815	-	rs41303421 Ref [C] Alt [A] Genotype1 [HOM_REF] Genotype2 [NO_CALL] Bases [C] Quals [10]
		// what's the missing data likelihood, is it right?
		char [] bases = {'C'};
		int [] qualities = {10};

	}


	@Test
	public void testImpossibleTwo () {
		char [] bases = {'A', 'A', 'T', 'T', 'T'};
		int [] qualities = {10,10,10,10, 10};
		VariantData vd = new VariantData(new Interval("chr1", 1, 1), 'A', 'T', GenotypeType.HOM_VAR, GenotypeType.HET, bases, qualities);

		int impossibleOne=vd.getNumImpossibleAllelesGenotypeOne();
		int impossibleTwo=vd.getNumImpossibleAllelesGenotypeTwo();
		Assert.assertEquals(impossibleOne, 2);
		Assert.assertEquals(impossibleTwo, 0);

	}
}
