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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.IntStream;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

/**
 * For a given set of bases and qualities, select the most common and 2nd most commonly seen alleles as ref/alt alleles.
 * Summarize the most common base quality in this pile-up by:
 * 1) Calculate the likelihood of the ref/ref, ref/alt, and alt/alt cases.
 * 2) Calculate the probability of the ref/ref class by taking the ref/ref likelihood and dividing it by the sum of all likelihoods (total probability).
 * 3) Convert this probability back to a Phread Score.
 * @author nemesh
 *
 */
public class SummarizeUMIBaseQualities {

	private final List<Byte> bases;
	private final List<Byte> qualities;
	private Byte mostCommonBase=null;
	// private Byte secondCommonBase=null;


	public SummarizeUMIBaseQualities (final List<Byte> bases, final List<Byte> qualities) {
		if (bases.size()!=qualities.size())
			throw new IllegalArgumentException("Number of bases [" + bases.size() +"] and qualities [" + qualities.size()+"] must be the same!");
		this.bases=bases;
		this.qualities=qualities;
	}

	public SummarizeUMIBaseQualities (final byte [] bases, final byte [] qualities) {
		this (convert(bases), convert(qualities));
	}

	/**
	 * Return the most commonly observed base in this pileup.
	 * Only looks at bases that occur, not at their quality.
	 * EXCEPTION - if there are exactly 2 observations and they disagree, 
	 * then the most common base is the one with the higher quality.
	 * @return
	 */
	public Byte getMostCommonBase() {
		// calculate the most common and second most common base in the same go
		// since you need to make the same object anyway.
		if (mostCommonBase!=null) return this.mostCommonBase;
		ObjectCounter<Byte> r = new ObjectCounter<>();
		for (Byte b : this.bases)
			r.increment(b);
		
		// this is a rare tie situation.  Select the highest quality base.
		if (getBasesTied(r)) {			
		    int maxIndex= IntStream.range(0, qualities.size())
		            .reduce((i, j) -> qualities.get(i) > qualities.get(j) ? i : j)
		            .getAsInt();
		    this.mostCommonBase=this.bases.get(maxIndex);		
		}  else {
			this.mostCommonBase=r.getMode();	
		}
		
		return this.mostCommonBase;
	}
	
	private boolean getBasesTied (ObjectCounter<Byte> r) {
		// if there's one key, it's not a tie.
		if (r.getKeys().size()==1) return false;
		// if there's more than 1 key, then see if all keys have the same number of counts.
		Set<Integer> counts = new HashSet<>(r.getCounts());
		if (counts.size()==1) return true; // all the same counts.
		return false;
	}
	
	/**
	 * Summarize the phread scores of all the bases piles in the given pileup.
	 * This calculates the mean quality score.
	 * For bases that disagree, this adds in the error rate.
	 * So, if you had AAT with qualities 30,30,10:
	 * Mean of : 0.999,0.999,0.1=0.6993333
	 * Thus an error rate of 1-0.6993333=0.3006667, which is a phread base score of ~ 5.
	 * We want to penalize bases with disagreements.
	 * Note: Phred scores of 0 are not allowed (can't have an error rate of 1), so if the mean is 0, it's set to 1.
	 * There's an issue here where the error rate will never be 1, but is quantized to 1 by the conversion to phred score.
	 * @return
	 */
	public int getSummarizedPhreadScoreByMeanWithErrors () {
		byte commonBase = getMostCommonBase();

		Mean mean = new Mean();

		for (int i=0; i<this.bases.size(); i++) {
			byte base = this.bases.get(i);
			byte qual = this.qualities.get(i);
			double prob = LikelihoodUtils.getInstance().phredScoreToErrorProbability(qual);

			if (base==commonBase)
				mean.increment(prob);
			else
				mean.increment(1-prob);
		}
		double meanErrorProbability = mean.getResult();
		int phred = LikelihoodUtils.getInstance().errorProbabilityToPhredScore(meanErrorProbability);
		// can't have a phread score of 0 for likeliylhood calculations.
		if (phred==0)
			phred=1;
		return phred;
	}

	/**
	 * Summarize the phread scores of all bases in the given pileup that are the most commonly observed base
	 * by the mode quality score.
	 * If there's more than 1 quality score that's the most common, the first is returned.
	 * @return
	 */
	public int getSummarizedPhreadScoreByMode () {
		byte commonBase = getMostCommonBase();
		ObjectCounter<Integer> baseQuals = new ObjectCounter<>();
		for (int i=0; i<this.bases.size(); i++) {
			byte base = this.bases.get(i);
			byte qual = this.qualities.get(i);
			if (base==commonBase)
				baseQuals.increment(Integer.valueOf(qual));
		}
		return baseQuals.getMode();
	}

	public int getSummarizedPhredScore () {
		return getSummarizedPhreadScoreByMeanWithErrors();
		//return getSummarizedPhreadScoreByMode();
	}

	private static List<Byte> convert (final byte [] d) {
		List<Byte> result = new ArrayList<>(d.length);
		for (byte b: d)
			result.add(b);
		return (result);
	}

}
