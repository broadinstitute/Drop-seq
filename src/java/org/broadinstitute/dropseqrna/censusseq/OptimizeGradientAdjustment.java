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

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.*;
import org.la4j.Vector;

/**
 * Given some current mixture and a gradient, how much should we apply the gradient (what scaling) to maximize the likelihood of the data?
 * @author nemesh
 *
 */
public class OptimizeGradientAdjustment implements UnivariateFunction {

	private double [] startingMixture;
	private double [] gradient;
	private OptimizeSampleRatiosLikelihoodFunctionCommonSNPs func;
	private double minimumMixture;
	private final int MAX_EVALUTATIONS=1000;

	public OptimizeGradientAdjustment (final double [] gradient, final double [] mixture, final OptimizeSampleRatiosLikelihoodFunctionCommonSNPs func, final double minimumMixture) {
		this.startingMixture=mixture;
		this.gradient=gradient;
		this.func=func;
		this.minimumMixture=minimumMixture;
	}

	@Override
	public double value(final double stepSize) {
		// update the mixture, calculate the likelihood.
		double [] newMix = updateMixtureWithGradient(gradient, startingMixture, stepSize);
		return func.value(newMix);
	}


	public UnivariatePointValuePair optimizeMixture (Double initialValue) {
		//UnivariateOptimizer minimizer = new BrentOptimizer(0.00001, 0.00001);
		UnivariateOptimizer minimizer = new BrentOptimizer(1e-8, 1e-12);
		// set the search interval.
		if (initialValue==null) initialValue=0.5;
		SearchInterval interval = new SearchInterval(0, 1, 0.5);
		// point is the mixture, value is the likelihood.
		UnivariatePointValuePair result = minimizer.optimize(new MaxEval(MAX_EVALUTATIONS), GoalType.MAXIMIZE, interval, new UnivariateObjectiveFunction(this));
		return (result);
	}



	double [] updateMixtureWithGradient(final double [] gradient, final double [] mixture, final double fudgeFactor) {
		Vector normlizedGradient = Vector.fromArray(gradient);
		// normalize to one, make distance moved smaller.
		// calculate the absolute sum.
		double absSum = 0;
		for (int i=0; i<normlizedGradient.length(); i++) {
			double absVal = Math.abs(normlizedGradient.get(i));
			absSum+=absVal;
		}
		normlizedGradient=normlizedGradient.divide(absSum);
		normlizedGradient=normlizedGradient.multiply(fudgeFactor);

		Vector finalMixture = Vector.fromArray(mixture);
		finalMixture=finalMixture.add(normlizedGradient);
		// normalize mixture to one.
		finalMixture= finalMixture.divide(finalMixture.sum());
		// check for out of bounds and correct.
		for (int i=0; i<finalMixture.length(); i++) {
			double val = finalMixture.get(i);
			if (val<minimumMixture)
				finalMixture.set(i, this.minimumMixture);
		}
		// normalize mixture to one, if values below minimum mixture were generated.
		finalMixture= finalMixture.divide(finalMixture.sum());
		// back to double []
		return finalMixture.toDenseVector().toArray();
	}

}
