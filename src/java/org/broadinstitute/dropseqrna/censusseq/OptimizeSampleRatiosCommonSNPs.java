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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;

import htsjdk.samtools.util.Log;

public class OptimizeSampleRatiosCommonSNPs {

	private static final Log log = Log.getInstance(OptimizeSampleRatiosCommonSNPs.class);

	private final OptimizeSampleRatiosLikelihoodFunctionCommonSNPs func;
	private final OptimizeSampleRatiosGradientFunction gradientFunc;
	private final CommonSNPsData data;
	private double convergenceMixtureThreshold=1e-4;
	private double minimumMixture=1e-12;
	private final boolean randomizedStarts;
	private final Random random;
	// if likelihood changes less than this amount, data converged.
	private final double minimumLikelihoodChange=0.1;

	public OptimizeSampleRatiosCommonSNPs(final CommonSNPsData data, final boolean scaleToDonorRepresentation, final int numThreads, final Random random) {
		this(data, scaleToDonorRepresentation, numThreads, random, false);
	}

	public OptimizeSampleRatiosCommonSNPs(final CommonSNPsData data, final boolean scaleToDonorRepresentation, final int numThreads, final Random random, final boolean randomizedStarts) {
		func = new OptimizeSampleRatiosLikelihoodFunctionCommonSNPs(data, numThreads);
		gradientFunc = new OptimizeSampleRatiosGradientFunction(data, scaleToDonorRepresentation, numThreads);
		this.data=data;
		this.randomizedStarts=randomizedStarts;
		this.random=random;
	}

	public OptimizeSampleRatiosCommonSNPsResult directIteration () {

		boolean converged=false;
		List<double [] > mixtureResults = new ArrayList<>();
		List<Double> likelihoodResults = new ArrayList<>();

		double[] startingMixture = getStartPoint(data.getSampleNames().size(), 1, this.randomizedStarts);
		mixtureResults.add(startingMixture);
		double startLikelihood=func.value(startingMixture);
		likelihoodResults.add (startLikelihood);
		log.info("STARTING CONDITIONS: Best likelihood [" + startLikelihood +"] mixture values "+ Arrays.toString(startingMixture));

		Double previousStepSize=null;
		// iterate.
		while (!converged) {
			int iteration=mixtureResults.size()-1;
			double [] mixture = mixtureResults.get(iteration);

			// attempt to find a step in the right direction.
			// calculate the gradient given the current mixture.
			double [] gradient=this.gradientFunc.value(mixture);
			// log.info("Gradient selected" + Arrays.toString(gradient));
			// what's the optimal step size along this gradient?
			OptimizeGradientAdjustment oga = new OptimizeGradientAdjustment(gradient, mixture, this.func, this.minimumMixture);
			UnivariatePointValuePair result = oga.optimizeMixture(previousStepSize);
			double stepSize = result.getPoint();
			double newLikelihood = result.getValue();
			//TODO uncomment this to test using the previous step size for the starting value of the next step.  Does this make optimization faster?
			previousStepSize=stepSize;

			// the new mixture given the optimized step size
			double [] newMixture = oga.updateMixtureWithGradient(gradient, mixture, stepSize);
			log.info("Best step size [" + stepSize +"] selected.  Best likelihood [" + newLikelihood +"] mixture values "+ Arrays.toString(newMixture));
			likelihoodResults.add(newLikelihood);
			mixtureResults.add(newMixture);

			// test if we should end.
			// Is the new likelihood better than a previous result?  Iteration moves from t0 to t1.
			iteration=mixtureResults.size()-1;

			double likelihood = likelihoodResults.get(iteration);
			double previousLikelihood = likelihoodResults.get(iteration-1);
			// if your likelihood is better than previously
			if (likelihood >= previousLikelihood) {
				double delta = likelihood-previousLikelihood;
				// and the likelihood hasn't changed much..
				if (delta<this.minimumLikelihoodChange) {
					log.info("Likelihood changed [" + delta +"] less than threshold ["+this.minimumLikelihoodChange +"]");
					converged=true;
				}
				else // or the mixtures haven't changed much...
					converged = mixturesConverged(mixtureResults, this.convergenceMixtureThreshold);
			}
			else {
				// you've overshot the best likelihood, stop.
				// maybe here you can explicitly drop donors to see if you can climb the hill more.
				log.info("Next application of gradient overshot best mixture, re-try or give up.");
				logGradientError(mixtureResults);
				break;
			}
		}
		// log how many iterations you went through.  This can be useful to note how "hard" it was to solve the mixture.
		log.info("Number of iterations before breaking out of convergence loop [", mixtureResults.size()-1+"]");

		// you've converged or the gradient is no longer helping. Get the best likelihood result and mixture.
		// useful if you overshot the best result (which can happen by 1 step.)
		double maxLikelihood =Double.NEGATIVE_INFINITY;
		int maxLikelihoodIndex=0;
		for (int i=0; i<likelihoodResults.size(); i++) {
			double val = likelihoodResults.get(i);
			if (val > maxLikelihood) {
				maxLikelihood=val;
				maxLikelihoodIndex=i;
			}
		}

		// now that you know which iteration did the best output the mixture results.
		double [] maximizedMixture = mixtureResults.get(maxLikelihoodIndex);

		// can this be reduced further by explicitly minimizing donors?


		List<String> donors =  func.getDonorNames();

		Map<String, Double> finalRatios = new HashMap<>();
		for (int i = 0; i < maximizedMixture.length; i++)
			finalRatios.put(donors.get(i), maximizedMixture[i]);

		// sort the likelihood results.
		Collections.sort(likelihoodResults);
		Collections.reverse(likelihoodResults);

		// TODO what if there was a polishing stage where each donor was reduced to the minimum mixture value and the likelihood calculated?
		// If there was any donor that was removed that made the data more likely, select the donor removal that improves likelihood score the most
		// Then iterate until there are no more donors to minimize.
		// TODO: when we have truly absent donors, are they the donors in the mixture that are fluctuating and thus not converging?
		// could we check all donors that were not removed for convergence separately?  If all the non-zero donors converge we're good?  (Can test if they had already converged?)
		// Would we need to run another round of optimization with those donors removed?  This might be more complicated if we have to remake the data object.
		OptimizeSampleRatiosCommonSNPsResult result = new OptimizeSampleRatiosCommonSNPsResult(finalRatios, converged, likelihoodResults.get(0), likelihoodResults.get(1), this.data.getTotalSNPAlleleCounts(), this.randomizedStarts);

		return result;
	}



	/**
	 * Have the mixtures only changed by a small amount?
	 * @param mixtureResults
	 * @return
	 */
	boolean mixturesConverged (final List<double [] > mixtureResults, final double threshold) {
		int index= mixtureResults.size()-1;
		double [] recentMixture = mixtureResults.get(index);
		double [] lastMixture = mixtureResults.get(index-1);
		boolean converged=true;
		for (int i=0; i<recentMixture.length; i++) {
			double rms=recentMixture[i];
			double lms=lastMixture[i];
			double delta = Math.abs(rms-lms);
			if (delta>threshold) {
				converged=false;
				break;
			}
		}
		if (converged)
			log.info("Mixtures converged! Delta less than threshold [" + threshold +"]");
			//log.info("New mixture" + Arrays.toString(recentMixture));
			//log.info("Old mixture" + Arrays.toString(lastMixture));
		return converged;
	}

	void logGradientError (final List<double [] > mixtureResults) {
		int index= mixtureResults.size()-1;
		double [] recentMixture = mixtureResults.get(index);
		double [] lastMixture = mixtureResults.get(index-1);
		String donorString = String.join(", ",func.getDonorNames());
		log.warn("Donor Names: [", donorString+"]");
		log.warn("Last good mixture result" + Arrays.toString(lastMixture));
		log.warn("Mis-step mixture result" + Arrays.toString(recentMixture));

	}


	double[] getStartPoint(final int n, final double value, final boolean randomizeStarts) {
		if (randomizeStarts) {
			double[] ds = new double[n];
			double total = 0;
			for (int i=0; i<ds.length;i++) {
				ds[i]=this.random.nextFloat();
				total+=ds[i];
			}
			for (int i=0; i<ds.length;i++)
				ds[i]=ds[i]/total;
			return (ds);
		} else {
			double frac = value / n;
			double[] ds = new double[n];
			Arrays.fill(ds, frac);
			return ds;
		}

	}

}
