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

import java.util.List;
import java.util.Map;

public class OptimizeSampleRatiosCommonSNPsResult {

	private final Map<String, Double> result;
	private final boolean converged;
	private final double bestLikelihood;
	private final double secondBestLikelihood;
	private final int numObservations;
	private final boolean randomizedStart;

	public OptimizeSampleRatiosCommonSNPsResult(final Map<String, Double> result, final boolean converged, final double bestLikelihood, final double secondBestLikelihood, final int numObservations, final boolean randomizedStart) {
		this.result=result;
		this.converged=converged;
		this.bestLikelihood = bestLikelihood;
		this.secondBestLikelihood=secondBestLikelihood;
		this.numObservations = numObservations;
		this.randomizedStart=randomizedStart;
	}

	public Map<String, Double> getResult() {
		return result;
	}
	
	public double [] getMixtureArray(List<String> donors) {
		double [] mixture=new double [donors.size()];
		for (int i=0; i<donors.size();i++)
			mixture[i]=result.get(donors.get(i));
		return mixture;						
	}

	public boolean isConverged() {
		return converged;
	}


	public double getBestLikelihood() {
		return this.bestLikelihood;
	}

	public int getNumObservations() {
		return numObservations;
	}

	/**
	 * Get the best likelihood normalized by the number of allelic observations.
	 * @return
	 */
	public double getNormalizedLikelihood () {
		return this.bestLikelihood/this.numObservations;
	}

	public double getSecondBestLikelihood() {
		return secondBestLikelihood;
	}

	public double getLikelihoodDelta() {
		return this.bestLikelihood-secondBestLikelihood;
	}

	public boolean isRandomizedStart() {
		return randomizedStart;
	}

}
