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
package org.broadinstitute.dropseqrna.utils.statistics;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;
import org.apache.commons.math3.stat.interval.ConfidenceInterval;
import org.apache.commons.math3.stat.interval.IntervalUtils;

/**
 * Calculates and stores the p-value and confidence interval for a set of observations.
 * Confidence intervals are calculated by the Wilson method
 * The current recommendations for the best CI method are the ones listed in the paper Interval Estimation for a Binomial Proportion by Brown, Cai and DasGupta in Statistical Science 2001, vol. 16, no. 2, pages 101-133. The authors examined several methods for calculating confidence intervals, and came to the following conclusion.
 * [W]e recommend the Wilson interval or the equal-tailed Jeffreys prior interval for small n and the interval suggested in Agresti and Coull for larger n.
 * We tend to have small samples sizes so Wilson it is.
 *
 * This is a wrapper class for convenience around the excellent org.apache.commons.math3 package.
 * See:
 * http://commons.apache.org/proper/commons-math/apidocs/org/apache/commons/math3/distribution/BinomialDistribution.html
 * https://commons.apache.org/proper/commons-math/apidocs/org/apache/commons/math3/stat/interval/IntervalUtils.html
 * https://commons.apache.org/proper/commons-math/apidocs/org/apache/commons/math3/stat/inference/AlternativeHypothesis.html
 * @author nemesh
 *
 */
public class BinomialStatistics {

	private final int numberOfTrials;
	private final int numberOfSuccesses;
	private final double ratio;
	private final ConfidenceInterval binomialConfidenceInterval;
	private final double binomialPvalue;
	// private final static DoubleRandomEngine rand = DoubleRandomEngine.makeDefault();

	/**
	 * Calculates binomial statistics for a data set. The assumed hypothesis is a two sided test.
	 * @param numberOfTrials The number of times you flipped a coin
	 * @param numberOfSuccesses The number of times you got a head (if you were testing for heads)
	 * @param probability The probability that the coin was weighted towards one side or the other.  Unbiased probability is 0.5.  Useful to set to a different value when
	 * alternativeHypothesis is set to a one sided test.
	 * @param confidenceLevel desired probability that the true probability of success falls within the returned interval.  This is typically 0.95 or 0.99.
	 */
	public BinomialStatistics(final int numberOfTrials, final int numberOfSuccesses, final double confidenceLevel, final double probability) {
		if (numberOfSuccesses>numberOfTrials)
			throw new IllegalArgumentException("Must have at least as many trials as successes.");
		this.numberOfTrials=numberOfTrials;
		this.numberOfSuccesses=numberOfSuccesses;
		this.ratio = (double) this.numberOfSuccesses / (double) this.numberOfTrials;

		this.binomialConfidenceInterval = IntervalUtils.getWilsonScoreInterval(numberOfTrials, numberOfSuccesses, confidenceLevel);

		// TODO: Which one of these is actually right?
		// truncated gives us the same values as the exact binomial test from R.
		this.binomialPvalue=getTwoSidedPvalueApacheTruncated(numberOfTrials, numberOfSuccesses, probability);
		//this.binomialPvalue=getTwoSidedPvalueApache(numberOfTrials, numberOfSuccesses, probability);

	}

	public static double getTwoSidedPvalueApacheTruncated(final int numTrials, final int numSuccesses, final double probability) {
		double pval = new BinomialTest().binomialTest(numTrials, numSuccesses, probability, AlternativeHypothesis.TWO_SIDED);
		if (pval>1) pval=1;
		return (pval);
	}

	/**
	public static double getTwoSidedPvalueColt(int numTrials, int numSuccesses, double probability) {
		double result =0;
		//  R code to replicate.
		//twoTailedUnbiasedTest<-function (numSuccesses, numTrials) {
		//	x=c(numSuccesses, numTrials-numSuccesses)
		//	dbinom(min(x), sum(x), 0.5) + ifelse( min(x)==0, 0, 2*sum(dbinom(0:(min(x)-1), sum(x), 0.5)))
		//


		// make numSuccesses less than 1/2 of the num successes.
		numSuccesses=Math.min(numSuccesses, numTrials-numSuccesses);

		//result=2*sum(dbinom(0:(min(x)), sum(x), 0.5)
		Binomial b = new Binomial(numTrials, probability, rand);

		for (int i=0; i<numSuccesses; i++) {
			result+=b.pdf(i);
		}

		result=result*2;
		result+=b.pdf(numSuccesses);

		return (result);
	}
	*/

	public static double getTwoSidedPvalue(final int numTrials, int numSuccesses, final double probability) {
		double result =0;

		/*  R code to replicate.
		twoTailedUnbiasedTest<-function (numSuccesses, numTrials) {
			x=c(numSuccesses, numTrials-numSuccesses)
			dbinom(min(x), sum(x), 0.5) + ifelse( min(x)==0, 0, 2*sum(dbinom(0:(min(x)-1), sum(x), 0.5)))
			#2*sum(dbinom(x=0:min(x[1], x[2]), size=x[1]+x[2], prob=0.5))
		}
		*/

		// make numSuccesses less than 1/2 of the num successes.
		numSuccesses=Math.min(numSuccesses, numTrials-numSuccesses);

		//result=2*sum(dbinom(0:(min(x)), sum(x), 0.5)
		BinomialDistribution bd = new BinomialDistribution(numTrials, probability);
		// double two = bd.cumulativeProbability(numberOfSuccesses);



		for (int i=0; i<numSuccesses; i++)
			result+=bd.probability(i);

		result=result*2;
		result+=bd.probability(numSuccesses);

		return (result);
	}


	public double getRatio() {
		return ratio;
	}

	public ConfidenceInterval getBinomialConfidenceInterval() {
		return binomialConfidenceInterval;
	}

	public double getBinomialPvalue() {
		return binomialPvalue;
	}

	@Override
	public String toString () {
		StringBuilder result = new StringBuilder();
		result.append("successes [" + this.numberOfSuccesses +"] trials [" + this.numberOfTrials + "] ratio [" + this.ratio +"] pval [" + this.binomialPvalue +"] ");
		result.append("confidence interval [" + this.binomialConfidenceInterval.getLowerBound() + "-" + this.binomialConfidenceInterval.getUpperBound() +"]");
		return result.toString();
	}


}
