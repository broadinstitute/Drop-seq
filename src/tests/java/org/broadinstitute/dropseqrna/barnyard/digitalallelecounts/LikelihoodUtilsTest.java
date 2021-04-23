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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.LikelihoodUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.GenotypeType;

/**
 * Note that the likelihoods have been changed from log to log 10.
 * @author nemesh
 *
 */
public class LikelihoodUtilsTest {

	
	@Test(enabled=true)
	public void testGetContaminationErrorRates () {
		// Easy first test.
		double [] r = LikelihoodUtils.getInstance().getContaminationErrorRates((byte) 20, 0.2, 0.2);
		Assert.assertEquals(r[0], 0.1684, 0.0001);
		Assert.assertEquals(r[1], 0.0496, 0.0001);
		
		// Test C=1, F=1 bound.
		r = LikelihoodUtils.getInstance().getContaminationErrorRates((byte) 20, 1, 1);
		Assert.assertEquals(r[0], 0.01, 0.0001);
		Assert.assertEquals(r[1], 0.99, 0.0001);
		
		
	}
	
	
	
	@Test(enabled=true)
	public void testMixedLikelihoodMultiRead () {
		GenotypeType [] g = {GenotypeType.HOM_REF, GenotypeType.HET, GenotypeType.HOM_VAR};
		List<GenotypeType> genotypes  = Arrays.asList(g);

		Double [] m = {new Double(2), new Double(1), new Double(1)};
		List<Double> mixture  = Arrays.asList(m);

		char refAllele ='A';
		char altAllele ='T';

		Byte [] b = {StringUtil.charToByte('A'), StringUtil.charToByte('A')};
		List<Byte> bases = Arrays.asList(b);
		Byte [] q = {new Byte ((byte)10), new Byte ((byte)10)};
		List<Byte> qualities =Arrays.asList(q);

		double result = LikelihoodUtils.getInstance().getLogLikelihoodMixedModel(refAllele, altAllele, genotypes, mixture, bases, qualities, null, null, null);
		Assert.assertEquals(result, Math.log10(0.36), 0.001);

	}

	@Test(enabled=true)
	public void testMixedLikelihood () {
		GenotypeType [] g = {GenotypeType.HOM_REF, GenotypeType.HET, GenotypeType.HOM_VAR};
		List<GenotypeType> genotypes  = Arrays.asList(g);

		Double [] m = {new Double(2), new Double(1), new Double(1)};
		List<Double> mixture  = Arrays.asList(m);

		char refAllele ='A';
		char altAllele ='T';

		List<Byte> bases = Collections.singletonList(StringUtil.charToByte('A'));
		List<Byte> qualities =Collections.singletonList(new Byte ((byte)10));

		double result = LikelihoodUtils.getInstance().getLogLikelihoodMixedModel(refAllele, altAllele, genotypes, mixture, bases, qualities, null, null, null);
		Assert.assertEquals(result, Math.log10(0.6), 0.001);
		
		double [] likes =  LikelihoodUtils.getInstance().getLikelihoodManyObservations ((byte) refAllele, (byte) altAllele, genotypes, bases.get(0), qualities.get(0), null, null, null, null, null);

		double result2=LikelihoodUtils.getInstance().getLikelihoodMixedModel(likes, mixture);
		Assert.assertEquals(Math.log10(result2), Math.log10(0.6), 0.001);
		
	}
	
	



	@Test(enabled=true)
	public void testLogLikelihood1 () {
		char [] b = {'A','A','A','A','A'};
		List<Byte> bases = convert (b);

		byte [] q = {10,10,10,10,10};
		List<Byte> qualities = convert (q);

		// probability if genotype is AA
		double result = LikelihoodUtils.getInstance().getLogLikelihood('A', 'A', bases, qualities, null, null, null, null, null);
		// > log10 (getRefLikelihood(5,0,0.9))
		// [1] -0.2287875
		double expected = -0.2287875;
		Assert.assertEquals(result,  expected, 0.000001);

		// probability if genotype is TT
		result = LikelihoodUtils.getInstance().getLogLikelihood('T', 'T', bases, qualities, null, null, null, null, null);
		// > log10 (getRefLikelihood(0,5,0.9))
		// [1] -5
		expected = -5.000002;
		Assert.assertEquals(result,  expected, 0.00001);

		// probability if genotype is AT
		result = LikelihoodUtils.getInstance().getLogLikelihood('A', 'T', bases, qualities, null, null, null, null, null);
		// log10(getHetLikelihood(5,0,0.9))
		// [1] -1.50515
		expected = -1.50515;
		Assert.assertEquals(result,  expected, 0.00001);
	}

	@Test(enabled=true)
	// true het genotype
	public void testLogLikelihood2 () {
		char [] b = {'A','A','A','T','T'};
		List<Byte> bases = convert (b);

		byte [] q = {20,20,20,20,20};
		List<Byte> qualities = convert (q);

		// probability if genotype is AA
		double result = LikelihoodUtils.getInstance().getLogLikelihood('A', 'A', bases, qualities, null, null, null, null, null);
		// log10 (getRefLikelihood(3,2,0.99))
		// [1] -4.013094
		double expected = -4.013094;
		Assert.assertEquals(result,  expected, 0.000001);

		// probability if genotype is TT
		result = LikelihoodUtils.getInstance().getLogLikelihood('T', 'T', bases, qualities, null, null, null, null, null);
		// log10 (getRefLikelihood(2,3,0.99))
		// [1] -6.00873
		expected = -6.008729;
		Assert.assertEquals(result,  expected, 0.00001);

		// probability if genotype is AT
		result = LikelihoodUtils.getInstance().getLogLikelihood('A', 'T', bases, qualities, null, null, null, null, null);
		// log10(getHetLikelihood(3,2,0.99))
		// [1] -1.50515
		expected = -1.50515;
		Assert.assertEquals(result,  expected, 0.00001);
	}

	@Test
	// Here, we copy the testLogLikelihood2 test, bump up the quality scores > 20, then max the probability at 20.
	public void testMaximumObservationProbability () {
		char [] b = {'A','A','A','T','T'};
		List<Byte> bases = convert (b);

		byte [] q = {30,30,30,30,30};
		List<Byte> qualities = convert (q);

		Double maxProb = 0.01;

		// probability if genotype is AA
		double result = LikelihoodUtils.getInstance().getLogLikelihood('A', 'A', bases, qualities, null, null, null, null, null);
		// log10 (getRefLikelihood(3,2,0.999))
		// [1] -6.001304
		double expected = -6.001304;
		Assert.assertEquals(result,  expected, 0.000001);
		// log10 (getRefLikelihood(3,2,0.99))
		// [1] -4.013094
		result = LikelihoodUtils.getInstance().getLogLikelihood('A', 'A', bases, qualities, null, maxProb, null, null, null);
		Assert.assertEquals(result,  -4.013094, 0.000001);

		// probability if genotype is TT
		result = LikelihoodUtils.getInstance().getLogLikelihood('T', 'T', bases, qualities, null, null, null, null, null);
		// log10 (getRefLikelihood(2,3,0.999))
		// [1] -9.000869
		expected = -9.000869;
		Assert.assertEquals(result,  expected, 0.00001);
		// log10 (getRefLikelihood(2,3,0.99))
		// [1] -6.00873
		result = LikelihoodUtils.getInstance().getLogLikelihood('T', 'T', bases, qualities, null, maxProb, null, null, null);
		Assert.assertEquals(result,  -6.00873, 0.000001);


		// probability if genotype is AT
		result = LikelihoodUtils.getInstance().getLogLikelihood('A', 'T', bases, qualities, null, null, null, null, null);
		// log10(getHetLikelihood(3,2,0.99))
		// [1] -1.50515
		expected = -1.50515;
		Assert.assertEquals(result,  expected, 0.00001);
		result = LikelihoodUtils.getInstance().getLogLikelihood('A', 'T', bases, qualities, null, maxProb, null, null, null);
		Assert.assertEquals(result,  expected, 0.00001); // doesn't matter what the error rate is for the het genotype state.

	}

	private List<Byte> convert (final char [] x) {
		List<Byte> r = new ArrayList<>();
		for (char c: x)
			r.add(new Byte((byte)c));
		return r;
	}

	private List<Byte> convert (final byte [] x) {
		List<Byte> r = new ArrayList<>();
		for (byte c: x)
			r.add(c);
		return r;
	}

	@Test(dataProvider = "PhreadToProb")
	public void testPhreadScoreToErrorProbability (final Byte phread, final Double prob) {
		double actual = LikelihoodUtils.getInstance().phredScoreToErrorProbability(phread);
		Assert.assertEquals(actual, prob, 0.0000001);
	}

	@Test
	public void testPhreadScoreToErrorProbability2 () {
		byte phread = 37;
		double prob=0.0001995262;
		double actual = LikelihoodUtils.getInstance().phredScoreToErrorProbability(phread);
		Assert.assertEquals(actual, prob, 0.00000001);
	}

	@Test(dataProvider = "PhreadToProb")
	public void testErrorProbabilityToPhreadScore (final Byte phread, final Double prob) {
		byte actual = LikelihoodUtils.getInstance().errorProbabilityToPhredScore(prob);
		Assert.assertEquals(actual, phread.byteValue());
	}

	@Test()
	public void testErrorProbabilityToPhreadScore2 () {
		double prob=0.0001995262;
		byte phread = 37;
		byte actual = LikelihoodUtils.getInstance().errorProbabilityToPhredScore(prob);
		Assert.assertEquals(actual, phread);
	}

	@Test
	public void testGetOneMinusPvalueFromLog10Likelihood () {
		LikelihoodUtils u = LikelihoodUtils.getInstance();
		byte [] phreds = new byte []  {60,50,40,35,30,20,10,5};
		double [] likes = new double [phreds.length];
		for (int i=0; i<phreds.length; i++ ) {
			likes[i]=  u.phredScoreToErrorProbability(phreds[i]);
			likes[i] = Math.log10(likes[i]);
		}

		double result = LikelihoodUtils.getInstance().getOneMinusPvalueFromLog10Likelihood(likes);
		double expected =0.26055;
		Assert.assertEquals(result, expected, 0.00001);

	}

	@Test
	public void testGetPvalueFromLog10Likelihood () {
		LikelihoodUtils u = LikelihoodUtils.getInstance();
		byte [] phreds = new byte []  {60,50,40,35,30,20,10,5};
		double [] likes = new double [phreds.length];
		for (int i=0; i<phreds.length; i++ ) {
			likes[i]=  u.phredScoreToErrorProbability(phreds[i]);
			likes[i] = Math.log10(likes[i]);
		}

		double resultOneMinus = LikelihoodUtils.getInstance().getOneMinusPvalueFromLog10Likelihood(likes);
		double expected =1-resultOneMinus;		
		double result=LikelihoodUtils.getInstance().getPvalueFromLog10Likelihood(likes);
		Assert.assertEquals(result, expected, 0.00001);

	}


	@DataProvider(name = "PhreadToProb")
	//https://en.wikipedia.org/wiki/Phred_quality_score
	public Object[][] createData1() {
	 return new Object[][] {
	   { new Byte((byte) 8), new Double(0.1584893) },
	   { new Byte((byte) 10), new Double(0.1) },
	   { new Byte((byte) 13), new Double(0.05011872) },
	   { new Byte((byte) 20), new Double(0.01) },
	   { new Byte((byte) 27), new Double(0.001995262) },
	   { new Byte((byte) 30), new Double(0.001) },
	   { new Byte((byte) 32), new Double(0.0006309573) },
	   { new Byte((byte) 37), new Double(0.0001995262) },
	   { new Byte((byte) 40), new Double(0.0001) },
	   { new Byte((byte) 50), new Double(0.00001) },
	   { new Byte((byte) 60), new Double(0.000001) }
	 };
	}
}
