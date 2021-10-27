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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment;

import java.util.Arrays;
import java.util.List;

import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.BestSampleAssignmentForCell;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.CellCollectionSampleLikelihoodCollection;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.CellSampleLikelihoodCollection;
import org.testng.Assert;
import org.testng.annotations.Test;

public class CellSampleLikelihoodCollectionTest {

	/**
	 * 	cell    		HUES62  				HUES53  				HUES64  				NKX     				HUES63  				ratio   				bestSample      trueSample
	 * 	GGCAGCCCTGAA    -416.628316553126       -135.091213037018       -504.328491670277       -627.125194530751       -475.01064616708        281.537103516108        HUES53  		HUES53
		AGAAAACGTGAT    -376.254865941466       -116.052455634381       -427.784172867756       -491.369446503429       -440.911250409602       260.202410307086        HUES53  		HUES53
		GCAGAATTAGGA    -16.3792016231805       -21.2208741400199       -21.2208741400199       -5.42459727705036       -18.069397703209        10.9546043461302        NKX     		NKX

		Results from the R implementation.
	 */

	@Test
	public void testAggregation () {
		List<String> donors = Arrays.asList("S1", "S2");
		
		CellCollectionSampleLikelihoodCollection ccslc = new CellCollectionSampleLikelihoodCollection(donors);
		
		CellSampleLikelihoodCollection aggregator = ccslc.buildCellSampleLikelihoodCollection("FOO");		

		CellSampleLikelihoodCollection c1 = ccslc.buildCellSampleLikelihoodCollection("FOO");
		c1.incrementNumSNPs(5);
		c1.incrementNumUMIs(6);
		c1.addLikelihood("S1", -3);
		c1.addLikelihood("S2", -4);

		CellSampleLikelihoodCollection c2 = ccslc.buildCellSampleLikelihoodCollection("FOO");
		c1.incrementNumSNPs(2);
		c1.incrementNumUMIs(3);
		c1.addLikelihood("S1", -8);
		c1.addLikelihood("S2", -6);

		aggregator.add(c1);
		aggregator.add(c2);

		double l1 = aggregator.getLoglikelihood("S1");
		Assert.assertEquals(l1, -11, 0.01);
		l1 = aggregator.getLoglikelihood("S2");
		Assert.assertEquals(l1, -10, 0.01);
		Assert.assertEquals(aggregator.getNumSNPs(), 7);
		Assert.assertEquals(aggregator.getNumUMIs(), 9);


	}

	@Test
	public void test3 () {
		List<String> donors = Arrays.asList("HUES62", "HUES53", "HUES64", "NKX", "HUES63");
		
		CellCollectionSampleLikelihoodCollection ccslc = new CellCollectionSampleLikelihoodCollection(donors);
		CellSampleLikelihoodCollection likelihoods = ccslc.buildCellSampleLikelihoodCollection("GCAGAATTAGGA");
		
		likelihoods.addLikelihood("HUES62", -16.3792016231805);
		likelihoods.addLikelihood("HUES53", -21.2208741400199);
		likelihoods.addLikelihood("HUES64", -21.2208741400199);
		likelihoods.addLikelihood("NKX", -5.42459727705036);
		likelihoods.addLikelihood("HUES63", -18.069397703209);

		BestSampleAssignmentForCell sampleAssignment = likelihoods.getBestSampleAssignment();
		double loglike = sampleAssignment.getBestLoglikelihood();
		double llRatio = sampleAssignment.getLogLikelihoodRatio();
		Assert.assertEquals(llRatio, 10.9546043461302, 0.000001);
		Assert.assertEquals(loglike, -5.42459727705036, 0.00001);
		Assert.assertEquals(sampleAssignment.getSample(), "NKX");
	}

	@Test
	public void test2 () {
		List<String> donors = Arrays.asList("HUES62", "HUES53", "HUES64", "NKX", "HUES63");
		CellCollectionSampleLikelihoodCollection ccslc = new CellCollectionSampleLikelihoodCollection(donors);
		CellSampleLikelihoodCollection likelihoods = ccslc.buildCellSampleLikelihoodCollection("AGAAAACGTGAT");
				
		likelihoods.addLikelihood("HUES62", -376.254865941466);
		likelihoods.addLikelihood("HUES53", -116.052455634381);
		likelihoods.addLikelihood("HUES64", -427.784172867756);
		likelihoods.addLikelihood("NKX", -491.369446503429);
		likelihoods.addLikelihood("HUES63", -440.911250409602);

		BestSampleAssignmentForCell sampleAssignment = likelihoods.getBestSampleAssignment();
		double loglike = sampleAssignment.getBestLoglikelihood();
		double llRatio = sampleAssignment.getLogLikelihoodRatio();
		Assert.assertEquals(llRatio, 260.202410307086, 0.000001);
		Assert.assertEquals(loglike, -116.052455634381, 0.00001);
		Assert.assertEquals(sampleAssignment.getSample(), "HUES53");
	}

	@Test
	public void test1 () {
		List<String> donors = Arrays.asList("HUES62", "HUES53", "HUES64", "NKX", "HUES63");
		CellCollectionSampleLikelihoodCollection ccslc = new CellCollectionSampleLikelihoodCollection(donors);
		CellSampleLikelihoodCollection likelihoods = ccslc.buildCellSampleLikelihoodCollection("GGCAGCCCTGAA");
						
		likelihoods.addLikelihood("HUES62", -416.628316553126);
		likelihoods.addLikelihood("HUES53", -135.091213037018);
		likelihoods.addLikelihood("HUES64", -504.328491670277);
		likelihoods.addLikelihood("NKX", -627.125194530751);
		likelihoods.addLikelihood("HUES63", -475.01064616708);

		BestSampleAssignmentForCell sampleAssignment = likelihoods.getBestSampleAssignment();
		double loglike = sampleAssignment.getBestLoglikelihood();
		double llRatio = sampleAssignment.getLogLikelihoodRatio();
		Assert.assertEquals(llRatio, 281.537103516108, 0.000001);
		Assert.assertEquals(loglike, -135.091213037018, 0.00001);
		Assert.assertEquals(sampleAssignment.getSample(), "HUES53");
	}





}
