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

import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;
import htsjdk.samtools.util.Log;
import jdistlib.Binomial;
import org.apache.commons.math3.analysis.MultivariateFunction;

import java.util.*;

public class OptimizeSampleRatiosLikelihoodFunctionPrivateSNPs implements MultivariateFunction {

	private final List<SNPSampleRecord> data;
	private static final Log log = Log.getInstance(OptimizeSampleRatiosLikelihoodFunctionPrivateSNPs.class);

	// in cases where the passed in ratio is 0.
	private final double MIN_RATIO = 1e-8;
	// in cases where the pvalue results to 0.
	private final double MIN_P_VALUE = 1e-64;

	// 0 based index from the sample name to the sample ratio.
	private final BiMap<String, Integer> sampleIndexMap;

	public OptimizeSampleRatiosLikelihoodFunctionPrivateSNPs(final List<SNPSampleRecord> data) {
		this.data = data;
		// build a map from the sample name to the index.
		sampleIndexMap = getSampleIndexMap(data);
	}

	/**
	 * When calling optimize on this function, the result will be a list of
	 *
	 * @return
	 */
	public Map<Integer, String> getRatioToSampleNameMap() {
		return sampleIndexMap.inverse();
	}

	public int numDonors() {
		return sampleIndexMap.size();
	}

	static BiMap<String, Integer> getSampleIndexMap(final List<SNPSampleRecord> data) {
		// first get a unique list of samples.
		Set<String> samples = new HashSet<>();

		for (SNPSampleRecord d : data)
			samples.add(d.getDonorName());

		BiMap<String, Integer> sampleIndexMap = HashBiMap.create();

		List<String> sortedSamples = new ArrayList<>(samples);
		Collections.sort(sortedSamples);

		int index = 0;
		for (String s : sortedSamples) {
			sampleIndexMap.put(s, index);
			index++;
		}
		return sampleIndexMap;
	}

	@Override
	// for a given set of sample ratios, what's the likelihood of observing all
	// the data?
	public double value(final double[] ratios) {
		double result = 0;
		double[] normalizedRatios = normalizeRatiosToOne(ratios);

		for (SNPSampleRecord r : this.data)
			result += Math.log10(evaluateSNPProbability(r, normalizedRatios));
		log.info("score: " + result + " param values " + Arrays.toString(normalizedRatios));
		return result;
	}

	/**
	 * For a SNP record and a set of ratios of samples, get the probability of
	 * this SNP's data given the ratios.
	 *
	 * @param r
	 * @param ratios
	 * @return
	 */
	private double evaluateSNPProbability(final SNPSampleRecord r, final double[] ratios) {

		int sampleIndex = sampleIndexMap.get(r.getDonorName()); // 0 based for
																// the array.
		double probability = ratios[sampleIndex];
		return getBinomialProb(r, probability);
	}

	public double getBinomialProb(final SNPSampleRecord r, final double probability) {
		int refCount = r.getRefCount();
		int altCount = r.getAltCount();
		int numTrials = refCount + altCount;
		double p = probability;
		if (r.getGenotype().equals("het"))
			p = p * 0.5;
		if (p < MIN_RATIO)
			p = MIN_RATIO;

		// BinomialDistribution bd = new BinomialDistribution(numTrials, p);
		// double score = bd.probability(altCount);
		// this is around 50x-100x faster than the apache commons package
		double score2 = Binomial.density(altCount, numTrials, p, false);
		if (score2 < this.MIN_P_VALUE)
			score2 = this.MIN_P_VALUE;

		return (score2);

	}

	public static double[] normalizeRatiosToOne(final double[] ratios) {
		double sum = 0;
		for (double ratio : ratios)
			sum += ratio;

		double[] result = new double[ratios.length];

		for (int i = 0; i < ratios.length; i++)
			result[i] = ratios[i] / sum;
		return (result);
	}

}
