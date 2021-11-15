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

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.LikelihoodUtils;

public class AllPairedSampleAssignmentsForCell {

	private SamplePairAssignmentForCell best;
	private List<SamplePairAssignmentForCell> all;
	private final String cell;
	private boolean dataAdded=false;
	private double bestPairPvalue;
	private final boolean scaleLikelihoods;
	
	public AllPairedSampleAssignmentsForCell (final String cell, final boolean scaleLikelihoods) {
		this.cell=cell;
		this.scaleLikelihoods=scaleLikelihoods;
		all=new ArrayList<>();
	}

	public void add (final SamplePairAssignmentForCell s) {				
		all.add(s);
		dataAdded=true;									
	}
	
	/**
	 * This looks for the pair of donors that is most likely to have generated the data
	 * It searches through all donor pairs where the results are well defined - they have at least one informative UMI to differentiate the two
	 * In the case where no donor pairs have at least one informative UMI (genetically similar donors with very few observations),
	 * select the pair with the most data, and set the mixture to be 0.5.  It's impossible to say what the correct mixture is, as all 3
	 * results are equally valid.
	 */
	private void calculateBest() {
		// set the scaling on the likelihoods of all donor pairs before recalculating.
		if (scaleLikelihoods) setScaling();
		
		// recalculate
		// TODO: should I instead compare all doublet likelihoods when computing the confident model doublet probability?
		// if the number of doublet likelihoods < 2, then the pvalue should be 1.
		// double [] doubleLikelihoods = all.stream().filter(x -> x.getDoubletPvalue()>=0.9d).mapToDouble(x -> x.getScaledDoubletLikelihood()).toArray();
		double [] likelihoods = all.stream().mapToDouble(x -> x.getScaledDoubletLikelihood()).toArray();
		double bestPairPvalue = LikelihoodUtils.getInstance().getPvalueFromLog10Likelihood(likelihoods);
		
		double maxLike=-Double.MAX_VALUE;
		SamplePairAssignmentForCell best = null;
		
		for (SamplePairAssignmentForCell spac: all) {
			double bestLike = spac.getScaledBestLikelihood();
			// important - you need to have informative UMIs for this to be a good call
			if (bestLike>maxLike & spac.getNumInformativeUMIs()>0) {
				maxLike=bestLike;
				best=spac;
			}
		}
						
		this.best=best;
		this.bestPairPvalue=bestPairPvalue;
		dataAdded=false;
		
		// fall back on finding the best uninformative pair.
		if (best==null) calculateBestUninformative();
	}
	
	/**
	 * Sets a scaling metric where donor pairs with fewer UMIs than the max have an additional penalty to their likelihoods to make them more consistent.
	 */
	private void setScaling () {
		int maxUMIs = all.stream().mapToInt(x -> x.getNumUMIs()).max().getAsInt();
		for (SamplePairAssignmentForCell s: all) {
			double scaling = (double) maxUMIs / s.getNumUMIs();
			s.setScalingFactor(scaling);
		}			
	}
	
	/**
	 * Scale likelihoods 
	 */
	
	/**
	 * In the case where no donor pairs have at least one informative UMI (genetically similar donors with very few observations),
	 * select the pair with the most data, and set the mixture to be 0.5.  It's impossible to say what the correct mixture is, as all 3
	 * results are equally valid.
	 * This can also result in the "best" pair having no UMIs at all.
	 */
	private void calculateBestUninformative() {
		double maxSize=0;
		SamplePairAssignmentForCell best = null;
		for (SamplePairAssignmentForCell spac: all) {
			int size = spac.getNumUMIs();
			if (size>=maxSize) {
				maxSize=size;
				best=spac;
			}
		}
		
		// construct a new object and set the mixture to be 0.5
		SamplePairAssignmentForCell result = new SamplePairAssignmentForCell(best.getCellBarcode(), best.getSampleOne(), best.getSampleTwo(), 
				best.getSampleOneSingleLikelihood(), best.getSampleTwoSingleLikelihood(), best.getDoubletLikelihood(), -1.0d, best.getImpossibleAllelesSampleOne(), 
				best.getImpossibleAllelesSampleTwo(), best.getNumInformativeSNPs(), best.getNumSNPs(), best.getNumUMIs(), best.getNumInformativeUMIs(), 
				best.getNumInformativeHomozygousUMIsSampleOne(), best.getNumInformativeHomozygousUMIsSampleTwo());
		this.best=result;
		
	}
	
	
	
	public String getCell() {
		return cell;
	}

	public SamplePairAssignmentForCell getBestAssignment () {
		// if data has been added, recalculate best.
		if (this.dataAdded) calculateBest();			
		return this.best;
	}
	
	public double getBestPairPvalue () {
		// if data has been added, recalculate best.
		if (this.dataAdded) calculateBest();			
		return this.bestPairPvalue;		
	}

	public List<SamplePairAssignmentForCell> getOtherAssignments() {
		List<SamplePairAssignmentForCell> result = new ArrayList<>(this.all);
		result.remove(this.best);
		return result;
	}

	/**
	 * Returns the results for a pair of donors.  If one of the donors given isn't in this data set, returns null.
	 * @param donorOne
	 * @param donorTwo
	 * @return
	 */
	public SamplePairAssignmentForCell getAssignmentForDonorPair(final String donorOne, final String donorTwo) {
		for (SamplePairAssignmentForCell o: this.all)
			if (o.getSampleOne().equals(donorOne) && o.getSampleTwo().equals(donorTwo))
				return o;
		return null;
	}
	
	


}
