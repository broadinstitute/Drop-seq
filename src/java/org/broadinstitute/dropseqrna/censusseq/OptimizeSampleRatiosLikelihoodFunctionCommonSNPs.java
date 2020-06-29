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

import htsjdk.samtools.util.Log;
import org.apache.commons.math3.analysis.MultivariateFunction;

import java.util.List;
import java.util.function.IntToDoubleFunction;
import java.util.stream.IntStream;

public class OptimizeSampleRatiosLikelihoodFunctionCommonSNPs implements MultivariateFunction {

	private final CommonSNPsData data;
	private static final Log log = Log.getInstance(OptimizeSampleRatiosLikelihoodFunctionCommonSNPs.class);
    private int numThreads;
    private Integer MAXIMUM_PENALITY=null;

    public OptimizeSampleRatiosLikelihoodFunctionCommonSNPs(final CommonSNPsData data, final int numThreads) {
		this.data=data;
		this.numThreads=numThreads;
		System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", Integer.toString(this.numThreads));
	}

	public OptimizeSampleRatiosLikelihoodFunctionCommonSNPs(final CommonSNPsData data) {
		this(data, 1);
	}

	@Override
	public double value(final double[] ratios) {
		double[] normalizedRatios = normalizeRatiosToOne(ratios);
		int numVariants = this.data.getNumVariants();

		IntToDoubleFunction calculateOne = (index) -> {
			int [] refAltCounts = this.data.getRefAltCounts(index);
			double minorAlleleFreq = this.data.getWeighedAlleleFrequenciesOneSNP(normalizedRatios, index);
			if (Double.isNaN(minorAlleleFreq))
				log.warn("NaN MAF detected!");
			double result = evaluateSNPProbability(minorAlleleFreq, refAltCounts);
			return result;
		};


		IntStream is = IntStream.range(0, numVariants);
		if (this.numThreads > 1)
			is = is.parallel();

		double result=is.mapToDouble(calculateOne).sum();

		//TODO: comment this out as it's pretty verbose.
		// log.info("score: " + result + " param values " + Arrays.toString(normalizedRatios));
		return result;

	}

	public int numDonors() {
		return data.getSampleNames().size();
	}

	public List<String> getDonorNames () {
		return this.data.getSampleNames();
	}

	/**
	 * The likelihood is:
	 * #a log (Fa) + b log (1-Fa)
	 * Where a is the count of the ref alleles, b is the count of the alt allelles, Fa is the frequency of the ref allele adjusted for the current
	 * sample mixture ratios.
	 * @param refAlleleFreq the reference allele freq
	 * @param refAltCounts
	 * @return
	 */
	public double evaluateSNPProbability(final double minorAlleleFreq, final int [] refAltCounts) {
		double result = (refAltCounts[0] * Math.log10(1-minorAlleleFreq)) + (refAltCounts[1] * Math.log10(minorAlleleFreq));
		// if a maximum penalty is set, use that instead of the likelihood when the likelihood penalty is stronger.
		if (this.MAXIMUM_PENALITY!=null)
			result=Math.max(result, this.MAXIMUM_PENALITY);
		if (Double.isNaN(result))
			log.warn("SNP results in NaN probability.  MAF [" + minorAlleleFreq +"] refCounts [" + refAltCounts[0] +"] altCounts [" +  refAltCounts[1]+"]");
		return result;
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
