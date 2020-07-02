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

import java.util.Arrays;
import java.util.function.BinaryOperator;
import java.util.function.IntFunction;
import java.util.stream.IntStream;

import org.apache.commons.math3.analysis.MultivariateVectorFunction;

import htsjdk.samtools.util.Log;

public class OptimizeSampleRatiosGradientFunction implements MultivariateVectorFunction {

	private final CommonSNPsData data;
	private static final Log log = Log.getInstance(OptimizeSampleRatiosLikelihoodFunctionCommonSNPs.class);
	private final int numThreads;

	public OptimizeSampleRatiosGradientFunction(final CommonSNPsData data) {
		this(data, 1);
	}

	public OptimizeSampleRatiosGradientFunction(final CommonSNPsData data, final int numThreads) {
		this.data=data;
		this.numThreads=numThreads;
		System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", Integer.toString(this.numThreads));
	}

	@Override
	public double [] value (final double [] ratios) {
		double[] normalizedRatios = normalizeRatiosToOne(ratios);
		int numVariants = this.data.getNumVariants();

		IntFunction<double []> calculateOne = (index) -> {
			int [] refAltCounts = this.data.getRefAltCounts(index);
			double minorAlleleFreq = this.data.getWeighedAlleleFrequenciesOneSNP(normalizedRatios, index);
			int [] genotypeStates = this.data.getCountsAltAllele(index);
			double [] result = new double [ratios.length];
			for (int j=0; j<result.length; j++)
				result[j]=gradientOneSNPOneSample(refAltCounts, minorAlleleFreq, genotypeStates[j]);
			return result;
		};

		BinaryOperator<double []> addResults= (a,b) -> {
			double [] result = new double [a.length];
			for (int i=0; i<a.length; i++)
				result[i]=a[i]+b[i];
			return result;
		};

		//TODO: this doesn't work in parallel stream context.  Is there a better way to reduce where we aren't generating new arrays constantly?
		BinaryOperator<double []> addResults2= (a,b) -> {
			for (int i=0; i<a.length; i++)
				a[i]=a[i]+b[i];
			return a;
		};

		double [] result = new double [ratios.length];
		Arrays.fill(result, 0);

		IntStream is = IntStream.range(0, numVariants);
		if (this.numThreads > 1)
			is = is.parallel();
		result=is.mapToObj(calculateOne).reduce(result, addResults);

		// double [] resultOld=valueOld(ratios);
		return result;
	}

	/**
	 * adjustmentFactor<-function (genotypeDonorState, numReadsRef, numReadsAlt, freqRef, freqAlt, weights=list("AA"=1, "AB"=0.5, "BB"=0)) {
     * 		if (!genotypeDonorState %in% names(weights)) warning(paste("Genotype State ilegal", genotypeDonorState))
     * 		if (freqRef==0 || freqAlt==0) warning ("Can't calculate for monomorphic SNPs!")
     * 		gi=weights[[genotypeDonorState]]
     * 		#if t1 is positive, there are more reference reads than expected.  We should upweigh AA samples.
     * 		t1=(numReadsRef/freqRef) - (numReadsAlt/freqAlt)
     * 		t2=(gi - freqRef)
     * 		return (t1*t2)
	 * }
	 *
	 * Calculate the gradient for the current mixture / genotype states.
	 * Missing genotype state values (-1) result in a gradient of 0 for that SNP.
	 *
	 * @return
	 */
	double gradientOneSNPOneSample (final int [] refAltCounts, final double minorAlleleFreq, final int genotypeState) {
		if (genotypeState==-1) return 0;
		double gs = (double) genotypeState/2;  // go from counts of alternate allele to fraction alternate allele.
		gs=Math.abs(1-gs);  //"FLIP" the 0 and 1 scores so back to the reference genotype is 1.
		double t1 = (refAltCounts[0]/(1-minorAlleleFreq)) - (refAltCounts[1]/(minorAlleleFreq));
		double t2 = (gs - (1-minorAlleleFreq));
		return t1*t2;
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
