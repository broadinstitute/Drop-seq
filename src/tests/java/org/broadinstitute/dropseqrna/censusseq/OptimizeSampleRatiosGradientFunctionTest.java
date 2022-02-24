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
package org.broadinstitute.dropseqrna.censusseq;

import org.broadinstitute.dropseqrna.censusseq.CommonSNPsData;
import org.junit.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;

public class OptimizeSampleRatiosGradientFunctionTest {
	/*
	 * $readCounts
  	refCount altCount
        9        1
        8        2
        6        4

	$genotypeStates
     		sample1 sample2 sample3 sample4
	[1,]     0.5     0.0     0.0     0.0
	[2,]     0.5     0.5     0.0     0.5
	[3,]     0.5     1.0     0.5     0.0

	sampleMixture 0.25 0.25 0.25 0.25
	minorAllelefreqs 0.125 0.375 0.500
	Gradient:
 	sample1   sample2   sample3   sample4
	-1.790476 -2.647619  3.085714  1.352381


	*/

	//TODO: should add some tests for donor scaling factor.  IE: OptimizeSampleRatiosGradientFunction(d, true)
	//This is tested on multiple large data sets to determine efficacy.
	@Test
	public void testGradient () {

		String[] sampleNames = { "sample1", "sample2", "sample3", "sample4" };
		int[][] gs = { { 1, 0, 0, 0 }, { 1, 1, 0, 1 }, { 1, 2, 1, 0 } };
		int[][] alleleCounts = { { 9, 1 }, { 8, 2 }, { 6, 4 } };
		// [{2, 4, 6}, {4, 7}, {6}, {4}]
		CommonSNPsData d = new CommonSNPsData(Arrays.asList(sampleNames));

		for (int i = 0; i < gs.length; i++) {
			int[] genos = gs[i];
			int[] ac = alleleCounts[i];
			d.addSNP(sampleNames, genos, ac[0], ac[1]);
		}

		OptimizeSampleRatiosGradientFunction f = new OptimizeSampleRatiosGradientFunction(d, false);
		double [] sampleMixtureStart = {0.25, 0.25, 0.25, 0.25};
		double [] result = f.value(sampleMixtureStart);
		double [] expected ={-1.790476, -2.647619,  3.085714,  1.352381};
		for (int i=0; i<result.length; i++)
			Assert.assertEquals(expected[i], result[i], 0.01);

	}


	/*
	 * 	For equal mixture:
	 *       sample1    sample2    sample3    sample4
			-2.9333333 -0.9333333  4.8000000  1.0666667
	 *
	 */
	@Test
	public void testGradientMissingValues () {

		String[] sampleNames = { "sample1", "sample2", "sample3", "sample4" };
		int[][] gs = { { 1, 0, 0, -1 }, { 1, 1, 0, 1 }, { 1, 2, 1, 0 } };
		int[][] alleleCounts = { { 9, 1 }, { 8, 2 }, { 6, 4 } };
		// [{2, 4, 6}, {4, 7}, {6}, {4}]
		CommonSNPsData d = new CommonSNPsData(Arrays.asList(sampleNames));

		for (int i = 0; i < gs.length; i++) {
			int[] genos = gs[i];
			int[] ac = alleleCounts[i];
			d.addSNP(sampleNames, genos, ac[0], ac[1]);
		}

		OptimizeSampleRatiosGradientFunction f = new OptimizeSampleRatiosGradientFunction(d, false);
		double [] sampleMixtureStart = {0.25, 0.25, 0.25, 0.25};
		double [] result = f.value(sampleMixtureStart);
		double [] expected ={-2.9333333, -0.9333333, 4.8000000, 1.0666667};
		for (int i=0; i<result.length; i++)
			Assert.assertEquals(expected[i], result[i], 0.01);

	}

}
