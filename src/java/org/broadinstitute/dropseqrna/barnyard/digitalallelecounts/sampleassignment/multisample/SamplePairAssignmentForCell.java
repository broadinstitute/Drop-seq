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

import org.apache.commons.math3.stat.StatUtils;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.LikelihoodUtils;

public class SamplePairAssignmentForCell {

	private final String cellBarcode;
	private final String sampleOne;
	private final String sampleTwo;

	// all likelihoods should be in log10 space.
	private final double sampleOneSingleLikelihood;
	private final double sampleTwoSingleLikelihood;
	private final double doubletLikelihood;

	// if a particular mixture of two cells is forced, keep that here.
	// private Double forcedMixture=null;
	// private Double forcedMixtureLikelihood=null;

	private final double mixture;
	private final int impossibleAllelesSampleOne;
	private final int impossibleAllelesSampleTwo;
	private final int numSNPs;
	private final int numInformativeSNPs;
	private final int numUMIs;
	private final int numInformativeUMIs;
	private final int numInformativeHomozygousUMIsSampleOne;
	private final int numInformativeHomozygousUMIsSampleTwo;
	
	// Normalize the likelihoods across many pairs of models for the same cell.
	// Likelihoods for each model are scaled by this factor.
	// This is used to select the best pair of donors out of all possible pairs, but has no other impact.
	private double scalingFactor;

	public SamplePairAssignmentForCell (final String cellBarcode, final String sampleOne, final String sampleTwo,
			final double sampleOneSingleLikelihood, final double sampleTwoSingleLikelihood, final double doubletLikelihood, final double mixture) {
		this.cellBarcode=cellBarcode;
		this.sampleOne=sampleOne;
		this.sampleTwo=sampleTwo;
		this.sampleOneSingleLikelihood=sampleOneSingleLikelihood;
		this.sampleTwoSingleLikelihood=sampleTwoSingleLikelihood;
		this.doubletLikelihood=doubletLikelihood;
		this.mixture=mixture;
		this.impossibleAllelesSampleOne=0;
		this.impossibleAllelesSampleTwo=0;
		this.numSNPs=0;
		this.numInformativeSNPs=0;
		this.numUMIs=0;
		this.numInformativeUMIs=0;
		this.numInformativeHomozygousUMIsSampleOne=0;
		this.numInformativeHomozygousUMIsSampleTwo=0;
		this.scalingFactor=1;
	}

	public SamplePairAssignmentForCell (final String cellBarcode, final String sampleOne, final String sampleTwo,
			final double sampleOneSingleLikelihood, final double sampleTwoSingleLikelihood, final double doubletLikelihood, final double mixture,
			final int impossibleAllelesSampleOne, final int impossibleAllelesSampleTwo, final int numInformativeSNPs, final int numSNPs, final int numUMIs, final int numInformativeUMIs,
			final int numInformativeHomozygousUMIsSampleOne, final int numInformativeHomozygousUMIsSampleTwo) {
		this.cellBarcode=cellBarcode;
		this.sampleOne=sampleOne;
		this.sampleTwo=sampleTwo;
		this.sampleOneSingleLikelihood=sampleOneSingleLikelihood;
		this.sampleTwoSingleLikelihood=sampleTwoSingleLikelihood;
		this.doubletLikelihood=doubletLikelihood;
		this.mixture=mixture;
		this.impossibleAllelesSampleOne=impossibleAllelesSampleOne;
		this.impossibleAllelesSampleTwo=impossibleAllelesSampleTwo;
		this.numInformativeSNPs=numInformativeSNPs;
		this.numSNPs=numSNPs;
		this.numUMIs=numUMIs;
		this.numInformativeUMIs=numInformativeUMIs;
		this.numInformativeHomozygousUMIsSampleOne=numInformativeHomozygousUMIsSampleOne;
		this.numInformativeHomozygousUMIsSampleTwo=numInformativeHomozygousUMIsSampleTwo;
		this.scalingFactor=1;
	}
	
	public double getSampleOneSingleLikelihood() {
		return sampleOneSingleLikelihood;
	}
	
	public double getScaledSampleOneSingleLikelihood () {
		return this.sampleOneSingleLikelihood*this.scalingFactor;
	}

	public double getSampleTwoSingleLikelihood() {
		return sampleTwoSingleLikelihood;
	}
	
	public double getScaledSampleTwoSingleLikelihood() {
		return this.sampleTwoSingleLikelihood*this.scalingFactor;
	}	

	public double getDoubletLikelihood() {
		return doubletLikelihood;
	}
	
	public double getScaledDoubletLikelihood() {
		return this.doubletLikelihood * this.scalingFactor;
	}

	public String getCellBarcode() {
		return cellBarcode;
	}

	public String getSampleOne() {
		return sampleOne;
	}

	public String getSampleTwo() {
		return sampleTwo;
	}

	public double getMixture() {
		return mixture;
	}

	public int getNumUMIs () {
		return this.numUMIs;
	}

	public int getImpossibleAllelesSampleOne() {
		return impossibleAllelesSampleOne;
	}

	public int getImpossibleAllelesSampleTwo() {
		return impossibleAllelesSampleTwo;
	}
	
	public int getNumInformativeHomozygousUMIsSampleOne() {
		return numInformativeHomozygousUMIsSampleOne;
	}

	public int getNumInformativeHomozygousUMIsSampleTwo() {
		return numInformativeHomozygousUMIsSampleTwo;
	}

	public int getNumSNPs() {
		return numSNPs;
	}

	public int getNumInformativeSNPs() {
		return numInformativeSNPs;
	}

	public double getDoubletLikelihoodRatio () {
		return this.doubletLikelihood - Math.max(this.sampleOneSingleLikelihood,  this.sampleTwoSingleLikelihood);
	}

	public double getBestLikelihood () {
		return Math.max(this.doubletLikelihood, Math.max(this.sampleOneSingleLikelihood, this.sampleTwoSingleLikelihood));
	}
	
	public double getScaledBestLikelihood () {
		return Math.max(getScaledDoubletLikelihood(), Math.max(getScaledSampleOneSingleLikelihood(), getScaledSampleTwoSingleLikelihood()));
	}
	
	/**
	 * Probability of a doublet is the likelihood of the doublet mixture divided by the sum of all likelihoods (doublet, singlet 1, singlet 2)
	 * @return The probability this cell is a doublet.
	 */
	public double getDoubletPvalue () {
		double [] allLikelihoods = {this.sampleOneSingleLikelihood, this.sampleTwoSingleLikelihood, this.doubletLikelihood};
		double maxValue = StatUtils.max(allLikelihoods);
		double totalLikelihood=0;
		
		for (int i=0; i<allLikelihoods.length; i++) {
			double d = allLikelihoods[i];
			d=d-maxValue;
			d=Math.pow(10, d);
			totalLikelihood+=d;
		}

		double mixtureLike = Math.pow(10, doubletLikelihood-maxValue);
		double result = mixtureLike/totalLikelihood;
		if (result < Double.MIN_VALUE) result = Double.MIN_VALUE;

		return result;
	}
	
	
	public String getCombinedDonorName() {
		return this.sampleOne+":"+this.sampleTwo;
	}

	public String getBestSample () {
		if (this.doubletLikelihood> Math.max(this.sampleOneSingleLikelihood, this.sampleTwoSingleLikelihood))
			return getCombinedDonorName();
		else if (this.sampleOneSingleLikelihood>this.sampleTwoSingleLikelihood)
			return this.sampleOne;
		else
			return this.sampleTwo;
	}
	
	public double getScalingFactor() {
		return scalingFactor;
	}

	public void setScalingFactor(double scalingFactor) {
		this.scalingFactor = scalingFactor;
	}

	

	@Override
	public String toString () {
		return this.cellBarcode +" Mixture [" + this.mixture+"] like [" + this.doubletLikelihood+"] S1 [" + this.sampleOne +"="+this.sampleOneSingleLikelihood+"] S2 ["
	+ this.sampleTwo+"="+this.sampleTwoSingleLikelihood+"]" + " impos alleles S1 ["+ this.impossibleAllelesSampleOne +"] impos alleles S2 [" + this.impossibleAllelesSampleTwo+"] "
	+ "num SNPs ["+this.numSNPs+"] num informative SNPs [" + this.numInformativeSNPs +"] num UMIs [" + this.numUMIs+"] num informative UMIs [" + this.numInformativeUMIs
	+ "] doublet likelihood ratio test [" + getDoubletLikelihoodRatio() +"] " + "bestLikelihood [" + this.getBestLikelihood() + "] best sample [" + this.getBestSample() +"]";
	}

	public int getNumInformativeUMIs() {
		return numInformativeUMIs;
	}
	
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((this.cellBarcode == null) ? 0 : cellBarcode.hashCode());
		long temp;
		temp = Double.doubleToLongBits(doubletLikelihood);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		result = prime * result + impossibleAllelesSampleOne;
		result = prime * result + impossibleAllelesSampleTwo;
		temp = Double.doubleToLongBits(mixture);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		result = prime * result + numInformativeSNPs;
		result = prime * result + numInformativeUMIs;
		result = prime * result + numSNPs;
		result = prime * result + numUMIs;
		result = prime * result + ((sampleOne == null) ? 0 : sampleOne.hashCode());
		temp = Double.doubleToLongBits(sampleOneSingleLikelihood);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		result = prime * result + ((sampleTwo == null) ? 0 : sampleTwo.hashCode());
		temp = Double.doubleToLongBits(sampleTwoSingleLikelihood);
		result = prime * result + (int) (temp ^ (temp >>> 32));
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		SamplePairAssignmentForCell other = (SamplePairAssignmentForCell) obj;
		if (cellBarcode == null) {
			if (other.cellBarcode != null)
				return false;
		} else if (!cellBarcode.equals(other.cellBarcode))
			return false;
		if (Double.doubleToLongBits(doubletLikelihood) != Double.doubleToLongBits(other.doubletLikelihood))
			return false;
		if (impossibleAllelesSampleOne != other.impossibleAllelesSampleOne)
			return false;
		if (impossibleAllelesSampleTwo != other.impossibleAllelesSampleTwo)
			return false;
		if (Double.doubleToLongBits(mixture) != Double.doubleToLongBits(other.mixture))
			return false;
		if (numInformativeSNPs != other.numInformativeSNPs)
			return false;
		if (numInformativeUMIs != other.numInformativeUMIs)
			return false;
		if (numSNPs != other.numSNPs)
			return false;
		if (numUMIs != other.numUMIs)
			return false;
		if (sampleOne == null) {
			if (other.sampleOne != null)
				return false;
		} else if (!sampleOne.equals(other.sampleOne))
			return false;
		if (Double.doubleToLongBits(sampleOneSingleLikelihood) != Double
				.doubleToLongBits(other.sampleOneSingleLikelihood))
			return false;
		if (sampleTwo == null) {
			if (other.sampleTwo != null)
				return false;
		} else if (!sampleTwo.equals(other.sampleTwo))
			return false;
		if (Double.doubleToLongBits(sampleTwoSingleLikelihood) != Double
				.doubleToLongBits(other.sampleTwoSingleLikelihood))
			return false;
		return true;
	}

	
	/*
	public Double getForcedMixture() {
		return forcedMixture;
	}

	public Double getForcedMixtureLikelihood() {
		return forcedMixtureLikelihood;
	}
	*/

	/**
	 * Compares the doublet likelihood to the two single likelihoods.
	 * If the doublet likelihood is equal to the sampleOneSingleLikelihood, then
	 * @param originalMixture
	 * @param doubletLikelihood
	 * @param sampleOneSingleLikelihood
	 * @param sampleTwoSingleLikelihood
	 */
	/*
	static double getMixture (final double originalMixture, final double doubletLikelihood, final double sampleOneSingleLikelihood, final double sampleTwoSingleLikelihood) {
		double epsilon = Math.min(doubletLikelihood,sampleOneSingleLikelihood)*0.01;
		int cmp = Precision.compareTo(doubletLikelihood, sampleOneSingleLikelihood, epsilon);
		// the first sample likelihood explains the data as well as the doublet likelihood.  The mixture is 1.
		if (cmp==0) return 1;

		epsilon = Math.min(doubletLikelihood,sampleTwoSingleLikelihood)*0.01;
		cmp = Precision.compareTo(doubletLikelihood, sampleTwoSingleLikelihood, epsilon);
		// the second sample likelihood explains the data as well as the doublet likelihood.  The mixture is 0.
		if (cmp==0) return 0;

		return originalMixture;

	}
	*/


}
