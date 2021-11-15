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

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.*;

import java.util.List;

public class VariantDataCollection implements UnivariateFunction {

	private final List<VariantData> vd;
	private final String sampleOne;
	private final String sampleTwo;
	private final String cell;

	public VariantDataCollection(final List<VariantData> vd, final String sampleOne, final String sampleTwo, final String cell) {
		this.vd=vd;
		this.sampleOne=sampleOne;
		this.sampleTwo=sampleTwo;
		this.cell=cell;
	}

	@Override
	public double value(final double mixture) {		
		double result=0;		
		for (VariantData vd: this.vd) {
			double like = vd.getLogLikelihood(mixture);
			result+=like;
		}		
		// double result = this.vd.parallelStream().mapToDouble(x -> x.getLogLikelihood(mixture)).sum();
		return result;
	}
		
	/**
	 * Generate a sample assignment for this cell, and add a likelihood value for the given mixture.
	 * @param forcedDonorMixture A mixture to set, can be null.
	 * @return
	 */
	public SamplePairAssignmentForCell optimizeMixture (final Double forcedDonorMixture) {
		UnivariateOptimizer minimizer = new BrentOptimizer(0.00001, 0.00001);

		// set the search interval.
		SearchInterval interval = null;
		if (forcedDonorMixture==null)
			interval= new SearchInterval(0, 1, 0.5);
		else
			interval=new SearchInterval(1-forcedDonorMixture, forcedDonorMixture, 0.5);

		// point is the mixture, value is the likelihood.
		UnivariatePointValuePair mixtureResult = minimizer.optimize(new MaxEval(100), GoalType.MAXIMIZE, interval, new UnivariateObjectiveFunction(this));
		double likelihoodSampleOne = value(1);
		double likelihoodSampleTwo=value(0);

		int numInformativeSNPs=getNumInformativeSNPs();
		int numUMIs=getNumUMIs();
		int numInformativeUMIs = getNumInformativeUMIs();
		
		double doubletLikelihood=mixtureResult.getValue();
		double mixture = mixtureResult.getPoint();
		// because the optimizer doesn't search the exact bounds of 0 and 1, need to correct the results to include those bounds.
		// if the first single sample explains the data better than the doublet, the mixture is 1.
		// don't do this if you're forcing a mixture interval to search around.
		if (likelihoodSampleOne>=doubletLikelihood && forcedDonorMixture==null) {
			doubletLikelihood=likelihoodSampleOne;
			mixture=1;
		}
		// if the second single sample explains the data better than the doublet, the mixture is 0.
		if (likelihoodSampleTwo>=doubletLikelihood && forcedDonorMixture==null) {
			doubletLikelihood=likelihoodSampleTwo;
			mixture=0;
		}

		SamplePairAssignmentForCell result = new SamplePairAssignmentForCell(this.cell, sampleOne, sampleTwo, likelihoodSampleOne, likelihoodSampleTwo, doubletLikelihood,
				mixture, getNumImpossibleAlleles(0), getNumImpossibleAlleles(1), numInformativeSNPs, vd.size(), numUMIs, numInformativeUMIs, 
				getNumInformativeHomozygousUMIs(0), getNumInformativeHomozygousUMIs(1));

		return result;
	}

	public SamplePairAssignmentForCell optimizeMixture () {
		return optimizeMixture(null);
	}

	private int getNumImpossibleAlleles(final int sampleIndex) {
		int count=0;
		for (VariantData v: this.vd)
			count+=v.getNumImpossibleAlleles(sampleIndex);
		return count;
	}

	private int getNumInformativeSNPs() {
		int count=0;
		for (VariantData v: this.vd)
			if(v.isInformative())
				count++;
		return count;
	}

	private int getNumInformativeUMIs() {
		int count=0;
		for (VariantData v: this.vd)
			if(v.isInformative())
				count+=v.getNumUMIs();
		return count;
	}
	
	private int getNumInformativeHomozygousUMIs(final int sampleIndex) {
		int count=0;
		for (VariantData v: this.vd)
			if(v.isInformativeHomozygous(sampleIndex))
				count+=v.getNumUMIs();
		return count;
	}
	

	private int getNumUMIs() {
		int count=0;
		for (VariantData v: this.vd)
			count+=v.getNumUMIs();
		return count;
	}

	public List<VariantData> getVariantData() {
		return this.vd;
	}

}
