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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.apache.commons.collections4.CollectionUtils;

/**
 * Given a cell's worth of data, a sample to test as donor 1, and a list of samples to test as other donors,
 * Test each sample pair and get a likelihood for the optimial mixture vs each individual in the pair.
 * @author nemesh
 *
 */
public class FindOptimalDonorMixture {

	private final VariantDataFactory factory;
	
	public FindOptimalDonorMixture(final VariantDataFactory f) {
		this.factory=f;
	}

	public AllPairedSampleAssignmentsForCell findBestDonorPair (final String sampleOne, final List<String> allSamples, final boolean scaleLikelihoods) {
		return findBestDonorPair(sampleOne, allSamples, null, scaleLikelihoods);
	}

	
	public AllPairedSampleAssignmentsForCell findBestDonorPair (final String sampleOne, final List<String> allSamples, final Double forcedMixture, final boolean scaleLikelihoods) {
		Collection<String> otherSamples = getNonPrimarySamples(sampleOne, allSamples);
		AllPairedSampleAssignmentsForCell result = new AllPairedSampleAssignmentsForCell(factory.getCell(), scaleLikelihoods);
		for (String sampleTwo: otherSamples) {
			VariantDataCollection vdc= factory.getVariantData(sampleOne, sampleTwo);
			SamplePairAssignmentForCell pair =  vdc.optimizeMixture(forcedMixture);
			result.add(pair);
		}		
		return result;
	}	
	
	 /* This is damn slow. */ 
	/* 
	public AllPairedSampleAssignmentsForCell findBestDonorPair (final String sampleOne, final List<String> allSamples, final Double forcedMixture) {
		Collection<String> otherSamples = getNonPrimarySamples(sampleOne, allSamples);
		AllPairedSampleAssignmentsForCell result = new AllPairedSampleAssignmentsForCell(factory.getCell());
		otherSamples.stream().forEach(x -> result.add(computePair(sampleOne, x, forcedMixture)));
		return result;
	}
	
	private SamplePairAssignmentForCell computePair (String sampleOne, String sampleTwo, final Double forcedMixture) {
		VariantDataCollection vdc= this.factory.getVariantData(sampleOne, sampleTwo);
		SamplePairAssignmentForCell pair =  vdc.optimizeMixture(forcedMixture);
		return (pair);
	}
	*/
	
	public static List<String> getNonPrimarySamples (final String sampleOne, final List<String> allSamples) {
		return new ArrayList<>(CollectionUtils.subtract(allSamples, Arrays.asList(sampleOne)));
	}
	
	
	

}
