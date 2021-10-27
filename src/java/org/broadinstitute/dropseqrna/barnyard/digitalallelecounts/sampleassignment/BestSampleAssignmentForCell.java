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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment;

/**
 * Simple container to hold the best sample assignment for a cell out of all possible assignments.
 * @author nemesh
 *
 */
public class BestSampleAssignmentForCell {

	private final String cellBarcode;
	private final String sample;
	private final double bestLoglikelihood;
	private final double secondBestLogLikelihood;
	private final double oneMinusPvalue;

	public BestSampleAssignmentForCell (final String cellBarcode, final String sample, final double bestLogLikelihood, final double secondBestLogLikelihood, final double oneMinusPvalue) {
		this.cellBarcode=cellBarcode;
		this.sample=sample;
		this.bestLoglikelihood=bestLogLikelihood;
		this.secondBestLogLikelihood=secondBestLogLikelihood;
		this.oneMinusPvalue=oneMinusPvalue;
	}

	public String getSample () {
		return this.sample;
	}

	public String getCellBarcode() {
		return cellBarcode;
	}

	/**
	 * Returns the best-secondBest loglikelihood
	 * @return
	 */
	public double getLogLikelihoodRatio () {
		return this.bestLoglikelihood-this.secondBestLogLikelihood;
	}

	public double getBestLoglikelihood () {
		return this.bestLoglikelihood;
	}

	public double getOneMinusPvalue() {
		return oneMinusPvalue;
	}


}
