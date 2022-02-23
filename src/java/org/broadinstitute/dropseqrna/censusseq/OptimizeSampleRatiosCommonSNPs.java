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
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.math3.optim.univariate.UnivariatePointValuePair;

import htsjdk.samtools.util.Log;

public class OptimizeSampleRatiosCommonSNPs {

	private static final Log log = Log.getInstance(OptimizeSampleRatiosCommonSNPs.class);

	private final OptimizeSampleRatiosLikelihoodFunctionCommonSNPs func;
	private final OptimizeSampleRatiosGradientFunction gradientFunc;
	private final CommonSNPsData data;
	private double convergenceMixtureThreshold=1e-4;
	private double minimumMixture=1e-5;
	// private double removalMixture=1e-8;
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

	private OptimizeSampleRatiosCommonSNPsResult directIteration (double[] startingMixture) {
		boolean converged=false;
		List<double [] > mixtureResults = new ArrayList<>();
		List<Double> likelihoodResults = new ArrayList<>();

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
		Map<String, Double> finalRatios=getRatioMap(maximizedMixture);
		
		// sort the likelihood results.
		Collections.sort(likelihoodResults);
		Collections.reverse(likelihoodResults);

		OptimizeSampleRatiosCommonSNPsResult result = new OptimizeSampleRatiosCommonSNPsResult(finalRatios, converged, likelihoodResults.get(0), likelihoodResults.get(1), this.data.getTotalSNPAlleleCounts(), this.randomizedStarts);

		return result;

	}
	
	private Map<String,Double> getRatioMap (double [] maximizedMixture) {
		List<String> donors =  func.getDonorNames();

		Map<String, Double> finalRatios = new HashMap<>();
		for (int i = 0; i < maximizedMixture.length; i++)
			finalRatios.put(donors.get(i), maximizedMixture[i]);
		return (finalRatios);
	}
	
	public OptimizeSampleRatiosCommonSNPsResult directIteration () {
		double[] startingMixture = getStartPoint(data.getSampleNames().size(), 1, this.randomizedStarts);
		OptimizeSampleRatiosCommonSNPsResult phaseOneResult = directIteration(startingMixture);
		return (phaseOneResult);
		
		// TODO: the likelihood calculation is undefined when the MAF is 0 or 1 - compromise by setting the donor to a very small fraction.
		// Need to solve that before this has even a chance of working.
		/*
		double [] roundOneMixture = phaseOneResult.getMixtureArray(this.func.getDonorNames());
		double [] zeroDonorOptimizedMixture =testForAbsentDonors(roundOneMixture);

		// if there were no changes made in the mixture, return the original optimization.
		if (Arrays.equals(roundOneMixture, zeroDonorOptimizedMixture)) {
			return (phaseOneResult);
		}
		
		Map<String, Double> finalRatios=getRatioMap(zeroDonorOptimizedMixture);
		
		OptimizeSampleRatiosCommonSNPsResult result = new OptimizeSampleRatiosCommonSNPsResult(finalRatios, phaseOneResult.isConverged(), phaseOneResult.getBestLikelihood(), 
				phaseOneResult.getSecondBestLikelihood(), this.data.getTotalSNPAlleleCounts(), this.randomizedStarts);
		
		return (result);
		*/				
	}
	
	/**
	 * After optimizing mixture, try to iteratively remove donors from the pool.
	 * Try to remove each donor and calculate the likelihood of the pool with the donor removed
	 * Select the donor that when removed best improves the total likelihood.
	 * If the likelihood does improve, set the donor proportion to 0, and test all non-zero donors
	 * Continue until there is no donor removal that will improve the pool's total likelihood.
	 *
	 * @param startingMixture The proportion of donors in the pool after a first pass optimization
	 * @return The proportion of donors after trying to remove donors from the pool 
	 */
	//TODO: fix likelihood at MAF=0 or MAF=1 to use this.  This is the "right" way to do this analysis.
	/*
	private double [] testForAbsentDonors (double [] startingMixture ) {
		
		
		double mixLikelihood=func.value(startingMixture);
		boolean searchActive=true;
		
		// make an explicit copy of the starting mixture.
		double [] mixture = startingMixture.clone();
				
		// while active, try to remove another donor each pass.
		while (searchActive) {
			int [] indexesToTest = getNonZeroRepDonorIndex(mixture);
			IntStream is = Arrays.stream(indexesToTest);
			// IntStream is = IntStream.range(0, startingMixture.length);
			//if (this.func.getNumThreads() > 1)
			//	is = is.parallel();
			log.info("Searching ["+indexesToTest.length+"] remaining donors (with representation <0.01) to remove donors who are not in pool.");
			List<Double> likelihoods=is.mapToObj(x-> getLikelihoodDonorRemoved(x,mixture)).collect(Collectors.toList());

			int index = likelihoods.indexOf(Collections.max(likelihoods));
			double bestLike = likelihoods.get(index);			
			// if you can't find a likelihood 
			if (bestLike<=mixLikelihood) {
				log.info("Unable to find any more donors to remove from pool.");
				break;
			}
				
			// improvement?
			if (bestLike>mixLikelihood) {
				// set the donor that improved things to mixture = 0;
				// since we're not testing all donors, we need to look up the correct original mixture index!
				int mixtureIndex=indexesToTest[index];
				mixture[mixtureIndex]=removalMixture;
				log.info("Donor at position [" + mixtureIndex+"] mixture set to ~0.  Previous pooled likelihood [" +mixLikelihood + "] new pooled likelihood [" + bestLike+"]");
				// update the best likelihood
				mixLikelihood=bestLike;
				
			}							
		}
		log.info("Finished removing donors from pool.");
		// the final optimized mixture.
		return mixture;
		
		
	}
	*/
	
	
	//TODO: this is a hackier way to do this, but the likelihood calculation is ill defined.
	/*
	private double [] testForAbsentDonors2 (double [] startingMixture ) {
		
		double mixLikelihood=func.value(startingMixture);
		log.info("Previous pooled likelihood [" +mixLikelihood + "]");
		
		// make an explicit copy of the starting mixture.
		double [] mixture = startingMixture.clone();
		
		int [] indexesToTest = getNonZeroRepDonorIndex(mixture);
		IntStream is = Arrays.stream(indexesToTest);
		
		
		List<Double> likelihoods=is.mapToObj(x-> getLikelihoodDonorRemoved(x,mixture)).collect(Collectors.toList());
		
		// get all indexes where donor removal would improve the results.
		for (int i=0; i<likelihoods.size(); i++) {
			double thisLike=likelihoods.get(i);
			
			if (likelihoods.get(i)>mixLikelihood) {
				log.info("Donor at position [" + i+"] mixture set to 0.  Pooled likelihood with this donor removed [" + thisLike+"]");
				mixture[i]=0;
			}
		}
		return mixture;
		
		
	}
	*/
	
	/**
	 * To save time, search for nonzero mixtures that are less than 0.01 mixture. 
	 * High mixture donors are unlikely to be truly absent from pool and reduce search space.
	 * @param mixture
	 * @return
	 */
	/*
	private int [] getNonZeroRepDonorIndex (double [] mixture) {
		List<Integer> result = new ArrayList<>();
		for (int i=0; i<mixture.length; i++) {
			if (mixture[i]>this.removalMixture & mixture[i]<0.01) result.add(i);
		}
		return result.stream().mapToInt(x->x).toArray();
	}
	*/
	
	/**
	 * For a given mixture vector, set the mixture at index to 0 and generate a likelihood.
	 * @param index The index of the donor to test for removal
	 * @param mixture The mixture vector
	 * @return the likelihood of the mixture with the donor removed.
	 */
	/*
	private Double getLikelihoodDonorRemoved (int index, double [] mixture) {
		// short circuit 
		if (mixture[index]==0) return Double.MIN_VALUE;
		double [] test = mixture.clone();
		test[index]=this.removalMixture;
		return func.value(test);						
	}
	/*


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
