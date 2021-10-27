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

import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.LikelihoodUtils;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.MultimapBuilder;

/**
 * For each cell, collect the likelihoods that the cell is taken from a particular sample
 * This collects the likelihoods across all samples, and stores the total likelihood across all SNPs
 * @author nemesh
 *
 */
public class CellSampleLikelihoodCollection {

	private final String cellBarcode;
	private int numUMIs;
	private int numSNPs;
	
	private final ImmutableMap<String, Integer> sampleIndexMap;
	private final double [] likelihoods;
	private double globalPenaltyScore;
	
	CellSampleLikelihoodCollection(final String cellBarcode, final int numSNPs, final int numUMIs, final double globalPenaltyScore, ImmutableMap<String, Integer> sampleIndexMap) {
		this.numSNPs=numSNPs;
		this.numUMIs=numUMIs;
		this.cellBarcode=cellBarcode;		
		this.sampleIndexMap=sampleIndexMap;
		this.likelihoods = new double [sampleIndexMap.size()];
		// explicitly set values to 0.
		Arrays.fill(this.likelihoods, 0);
		this.globalPenaltyScore=globalPenaltyScore;
	}
	
	void add (final CellSampleLikelihoodCollection other) {
		if (!other.getCellBarcode().equals(this.cellBarcode))
			throw new IllegalArgumentException("Adding a collection in with the wrong cell barcode.  Original ["+this.cellBarcode+"] new ["+ other.getCellBarcode()+"]");
		this.numSNPs+=other.getNumSNPs();
		this.numUMIs+=other.getNumUMIs();
		for (String sample: this.sampleIndexMap.keySet())
			addLikelihood(sample, other.getLoglikelihood(sample));
		this.globalPenaltyScore+=other.getGlobalPenaltyScore();

	}
	
	void addLikelihood(final String sample, final double logLikelihood) {
		int index = this.sampleIndexMap.get(sample);
		this.likelihoods[index]+=logLikelihood;
	}
	
	/**
	 * Increase the global likelihood penalty score by this log likelihood.
	 * @param logLikelihood
	 */
	void addGlobalPenalty(final double logLikelihood) {
		this.globalPenaltyScore+=logLikelihood;
	}
		
	/**
	 * Get the median likelihood across all donors.
	 * @return
	 */
	public double getMedianLikelihood () {
		Median median = new Median();
		return median.evaluate(likelihoods);
	}

	public String getCellBarcode () {
		return this.cellBarcode;
	}

	public void incrementNumSNPs (final int num) {
		this.numSNPs+=num;
	}

	public void incrementNumUMIs (final int num) {
		this.numUMIs+=num;
	}

	public int getNumUMIs() {
		return numUMIs;
	}

	public int getNumSNPs() {
		return numSNPs;
	}


	/**
	 * Returns the log likelihood gathered for this sample, or null if this sample isn't in the data set.
	 * @param sample
	 * @return
	 */
	public Double getLoglikelihood (final String sample) {	
		Integer index=this.sampleIndexMap.get(sample);
		if (index==null) return null;
		return this.likelihoods[index];
	}

	public boolean isSampleInCollection (final String sample) {
		return this.sampleIndexMap.containsKey(sample);
	}

	public Set<String> getSamples () {
		return this.sampleIndexMap.keySet();
	}

	
	/**
	 * Based on all data added, computes the best log likelihood and 2nd best likelihood sample assignments.
	 * @return An object containing the best sample, 1st and 2nd best log likelihoods.
	 */
	// 
	public BestSampleAssignmentForCell getBestSampleAssignment () {
		// track best/second best likelihoods.
		Double bestLogLike = Double.NEGATIVE_INFINITY;
		Double secondBestLogLike=Double.NEGATIVE_INFINITY;
		String bestSample=null;
				
		for (String sample: sampleIndexMap.keySet()) {			
			double val = this.getLoglikelihood(sample);
			if (val>bestLogLike) {
				// make the best the 2nd best.
				secondBestLogLike=bestLogLike;
				bestSample=sample;
				bestLogLike=val;								
			} else {
				if (val>secondBestLogLike) {
					secondBestLogLike=val;
				}
			}
		}
						
		double pval = LikelihoodUtils.getInstance().getOneMinusPvalueFromLog10Likelihood(this.likelihoods);
		BestSampleAssignmentForCell best = new BestSampleAssignmentForCell(this.cellBarcode, bestSample, bestLogLike, secondBestLogLike, pval);
		return best;
	}
	
	public List<String> getDonorsRankedByAssignmentLikelihood () {
		String[] sampleNames = sampleIndexMap.keySet().stream().toArray(String []::new);
		
		List<String> sortedNames = IntStream.range(0, likelihoods.length)
			    .mapToObj(i -> i)
			    .sorted(Comparator.comparing(i -> likelihoods[i], Comparator.reverseOrder()))
			    .map(i -> sampleNames[i])
			    .collect(Collectors.toList());
		
		return sortedNames;
		
		/*
		ListMultimap<Double, String> treeListMultimap = MultimapBuilder.treeKeys().arrayListValues().build();
		for (String donor: sampleIndexMap.keySet()) {
			int index = sampleIndexMap.get(donor);
			double value = this.likelihoods[index];
			treeListMultimap.put(value, donor);
		}
		*/
		// map sorted by key already?  Which direction?
		
		
		
		
		
	}

	public double getGlobalPenaltyScore() {
		return globalPenaltyScore;
	}
	
}
