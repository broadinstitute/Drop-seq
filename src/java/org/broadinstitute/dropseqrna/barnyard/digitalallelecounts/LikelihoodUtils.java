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

import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.stat.StatUtils;

import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.GenotypeType;

/**
 * Helper methods for calculating likelihoods and converting to/from Phred scores.
 * @author nemesh
 *
 */
public class LikelihoodUtils {


	private static LikelihoodUtils instance = null;
	// holds a list of phred to probability scores, ie: 10->0.9, 20->0.99, 30->0.999, etc.
	private double [] phredToProbability;

	private LikelihoodUtils() {
		phredToProbability=new double [128];
		Arrays.fill(phredToProbability, -1d);
	}

	public static synchronized LikelihoodUtils getInstance() {
        if (instance == null)
			instance = new LikelihoodUtils();

        return instance;
    }


	/**
	 * Given a set of likelihoods in log10, output 1- the probability of the most largest likelihood [ 1-p].
	 * @param allLikelihoods a collection of likelihoods
	 * @return 1 - (best like / sum (likes))
	 */
	public double getOneMinusPvalueFromLog10Likelihood (final double [] allLikelihoods) {
		// we clone the array so we don't change it.
		double [] likes=allLikelihoods.clone();
		Arrays.sort(likes);

		double maxValue = StatUtils.max(likes);

		double totalLikelihood=0;
		double allButBestLikelihood=0;

		for (int i=0; i<likes.length; i++) {
			double d = likes[i];
			d=d-maxValue;
			d=Math.pow(10, d);
			totalLikelihood+=d;
			if (i!=(likes.length-1))
				allButBestLikelihood+=d;
		}

		double result = allButBestLikelihood/totalLikelihood;
		if (result < Double.MIN_VALUE) result = Double.MIN_VALUE;
		return result;
	}
	
	/**
	 * Given a set of likelihoods in log10, output the probability of the most largest likelihood [p].
	 * @param allLikelihoods
	 * @return
	 */
	public double getPvalueFromLog10Likelihood (final double [] allLikelihoods) {
		//TODO: is it better to clone the array here, or should the class handing off the array do the cloning?
		double [] likes=allLikelihoods.clone();
		Arrays.sort(likes);
	
		double maxValue = StatUtils.max(likes);
		double totalLikelihood=0;
		
		for (int i=0; i<likes.length; i++) {
			double d = likes[i];
			d=d-maxValue;
			d=Math.pow(10, d);
			totalLikelihood+=d;
		}

		double result = 1/totalLikelihood;
		if (result < Double.MIN_VALUE) result = Double.MIN_VALUE;
		return result;
	}
	
	

	/**
	 * Calculate the likelihood for a pileup of bases and qualities using a mixture of different models.
	 * For example, one could calculate a mixture of two genotypes, a ref and a het, by adding those two models to the states list, and adding to 0.5 entries to the mixture list.
	 * For each base, the probability of each of the genotypes is calculated and multiplied by the mixture for that genotype.  The individual likelihoods are then summed,
	 * and the result is then divided by the sum of the mixtures.
	 * This is done per observation, then the log is taken and the results of all log-likelihoods are summed.
	 *
	 * @param refAllele the reference allele for the variant
	 * @param altAllele the alternate allele for the variant
	 * @param The list of genotype states for the observed genotypes for the variant (donors that correctly genotyped)
	 * @param mixture The list of mixture parameters that these genotypes are mixed by.  This list must be in the same order as the genotype states.
	 * @param bases The bases observed
	 * @param qualities The qualities of the bases observed.  This list must be in the same order as the bases observed.
	 * @return
	 */
	public double getLogLikelihoodMixedModel (final char refAllele, final char altAllele, final List<GenotypeType> genotypes,
			final List<Double> mixture, final List<Byte> bases, final List<Byte> qualities, final Double missingDataPenality, final Double genotypeProbability, final Double maximumObservationProbability) {

		byte ref= StringUtil.charToByte(refAllele);
		byte alt= StringUtil.charToByte(altAllele);
		double result = getLogLikelihoodMixedModel(ref, alt, genotypes, mixture, bases, qualities, missingDataPenality, genotypeProbability, maximumObservationProbability);
		return result;
	}

	public double getLogLikelihoodMixedModel (final byte refAllele, final byte altAllele, final List<GenotypeType> genotypes,
			final List<Double> mixture, final List<Byte> bases, final List<Byte> qualities, final Double missingDataPenality, final Double genotypeProbability, final Double maximumObservationProbability) {

		if (genotypes.size()!=mixture.size())
			throw new IllegalArgumentException("Genotype list and mixture list must be the same size.");

		double result=0;

		for (int i=0; i<bases.size(); i++) {
			double likelihood = getLikelihoodMixedModel(refAllele, altAllele, genotypes, mixture, bases.get(i), qualities.get(i), missingDataPenality, genotypeProbability, maximumObservationProbability);			
			result+=Math.log10(likelihood);						
		}
		return result;
	}
	
	public double getLogLikelihoodMixedModel (double [] [] likelihoods, final List<Double> mixture) {

		if (likelihoods[0].length!=mixture.size())
			throw new IllegalArgumentException("likelihood list and mixture list must be the same size.");

		double result=0;

		for (int i=0; i<likelihoods.length; i++) {
			double likelihood = getLikelihoodMixedModel(likelihoods[i], mixture);			
			result+=Math.log10(likelihood);						
		}
		return result;
	}


	/**
	 * Calculates the likelihood for a list of genotype states for a single observation.
	 * @param ref the reference allele for the variant
	 * @param alt the alternate allele for the variant
	 * @param genotypes The list of genotypes to calculate likelihoods for
	 * @param mixture A list of mixture coefficientss for the list of genotypes.
	 * @param base The observed base in the sequence pileup
	 * @param quality The quality of the observed base
	 * @param missingDataPenality a pre-computed likelihood to use instead of computing a likelihood for a genotype.
	 * @param genotypeProbability The likelihood of the genotype.  Can be null to ignore.
	 * @param maximumObservationProbability If set, this iss the maximum penalty that can be generated for a single observation.
	 * @return the likelihood of the base/quality, given a mixture of genotypes.
	 */
	private double getLikelihoodMixedModel(final byte ref, final byte alt, final List<GenotypeType> genotypes, final List<Double> mixture,
			final Byte base, final Byte quality, final Double missingDataPenality, final Double genotypeProbability, final Double maximumObservationProbability) {

		double result=0;
		double sumMixture=0;
		// for each genotype model, calculate the likelihood, multiply by the mixture for that model, then add to the result.
		for (int i=0; i<genotypes.size(); i++) {
			GenotypeType genotype = genotypes.get(i);
			if ((genotype==GenotypeType.NO_CALL || genotype==GenotypeType.UNAVAILABLE) && missingDataPenality==null)
				throw new IllegalArgumentException("If using NO_CALL or UNAVAILABLE genotypes, must set a missingDataPenality!");
			double mix = mixture.get(i);
			sumMixture+=mix;
			double likelihood = 0;
			// calculate the single likelihood.
			switch (genotype) {
			case HOM_REF:
				likelihood=getLikelihoodHomozygote(ref, ref, base, quality, null, maximumObservationProbability); break;
			case HOM_VAR:
				likelihood=getLikelihoodHomozygote(alt, alt, base, quality, null, maximumObservationProbability); break;
			case HET:
				likelihood=getLikelihoodHeterozygote(ref, alt, base, quality, null, maximumObservationProbability); break;
			case NO_CALL:
				likelihood=missingDataPenality; break;
			case UNAVAILABLE:
				likelihood=missingDataPenality; break;
			default:
			}
			// multiply the calculated model likelihood by the ratio
			likelihood=likelihood*mix;
			// multiply the result by the genotype probability if it is not null
			if (genotypeProbability!=null)
				likelihood*=genotypeProbability;
			result+=likelihood;
		}
		result=result/sumMixture;		
		return result;
	}
	
	/**
	 * Calculates the likelihood for a list of genotype states for a single UMI observation.
	 * Here, the likelihoods for the genotype states have been precomputed.  The number of likelihoods should be equal to the number of mixture coefficients.
	 * This is the equivilent of the weighted average of the likelihoods.
	 * @param likelihoods An array of doubles representing the likelihoods of the genotypes in this mixture
	 * @param mixture A list of mixture coefficients to weight the genotype likelihoods.
	 * @return The likelihood of the mixed model
	 * @see getLikelihoodManyObservations
	 */
	public double getLikelihoodMixedModel(double [] likelihoods, final List<Double> mixture) {
		double result=0;
		double sumMixture=0;
		for (int i=0; i<likelihoods.length; i++) {
			double mix = mixture.get(i);
			sumMixture+=mix;
			double likelihood = likelihoods[i]*mixture.get(i);
			result+=likelihood;
		}
		result=result/sumMixture;
		return result;
	}
	
	/**
	 * Get a list of likelihoods for a given base/quality for a list of genotypes.
	 * @param ref the reference allele for the variant
	 * @param alt the alternate allele for the variant
	 * @param genotypes The list of genotypes to calculate likelihoods for
	 * @param base The observed base in the sequence pileup
	 * @param quality The quality of the observed base
	 * @param missingDataPenality a pre-computed likelihood to use instead of computing a likelihood for a genotype.
	 * @param genotypeProbability The likelihood of the genotype.  Can be null to ignore.
	 * @param maximumObservationProbability If set, this iss the maximum penalty that can be generated for a single observation.
	 * @return An array of likelihoods in the same order as the submitted genotypes.
     */
	public double [] getLikelihoodManyObservations (final byte ref, final byte alt, final List<GenotypeType> genotypes, final Byte base, final Byte quality, 
			final Double missingDataPenality, final Double genotypeProbability, final Double maximumObservationProbability) {
		
		double [] result = new double [genotypes.size()];
		
		for (int i=0; i<genotypes.size(); i++) {
			GenotypeType genotype = genotypes.get(i);
			if ((genotype==GenotypeType.NO_CALL || genotype==GenotypeType.UNAVAILABLE) && missingDataPenality==null)
				throw new IllegalArgumentException("If using NO_CALL or UNAVAILABLE genotypes, must set a missingDataPenality!");
			double likelihood = 0;
			// calculate the single likelihood.
			switch (genotype) {
			case HOM_REF:
				likelihood=getLikelihoodHomozygote(ref, ref, base, quality, null, maximumObservationProbability); break;
			case HOM_VAR:
				likelihood=getLikelihoodHomozygote(alt, alt, base, quality, null, maximumObservationProbability); break;
			case HET:
				likelihood=getLikelihoodHeterozygote(ref, alt, base, quality, null, maximumObservationProbability); break;
			case NO_CALL:
				likelihood=missingDataPenality; break;
			case UNAVAILABLE:
				likelihood=missingDataPenality; break;
			default:
			}
			result[i]=likelihood;
		}				
		return result;
		
	}

	public double getLogLikelihood (final char refAllele, final char altAllele, final List<Byte> bases, final List<Byte> qualities, final Double genotypeQuality, final Double maximumObservationProbability) {
		byte ref= StringUtil.charToByte(refAllele);
		byte alt= StringUtil.charToByte(altAllele);
		return getLogLikelihood(ref, alt, bases, qualities, genotypeQuality, maximumObservationProbability);
	}

	public double getLogLikelihood (final byte refAllele, final byte altAllele, final List<Byte> bases, final List<Byte> qualities, final Double genotypeQuality, final Double maximumObservationProbability) {

		if (refAllele==altAllele)
			return getLogLikelihoodHomozygote(refAllele, altAllele, bases, qualities, genotypeQuality, maximumObservationProbability);
		else
			return getLogLikelihoodHeterozygote(refAllele, altAllele, bases, qualities, genotypeQuality, maximumObservationProbability);
	}

	public double getLikelihood (final byte refAllele, final byte altAllele, final Byte base, final Byte quality, final Double genotypeQuality, final Double maximumObservationProbability) {
		if (refAllele==altAllele)
			return getLikelihoodHomozygote(refAllele, altAllele, base, quality, genotypeQuality, maximumObservationProbability);
		else
			return getLikelihoodHeterozygote(refAllele, altAllele, base, quality, genotypeQuality, maximumObservationProbability);
	}


	/**
	 * The one or more observation version to get the homozygous log likelihood.
	 * Iteratively calculates the homozygous log likelihood, takes the log, and sums those logs.
	 * @param refAllele
	 * @param altAllele
	 * @param bases
	 * @param qualities
	 * @return
	 */
	private double getLogLikelihoodHomozygote (final byte refAllele, final byte altAllele, final List<Byte> bases, final List<Byte> qualities, final Double genotypeQuality, final Double maximumObservationProbability) {
		if (refAllele!=altAllele)
			throw new IllegalArgumentException("For homozygous likelihood, ref allele [" + refAllele +"] and alt allele [" + altAllele+ "] must match!");

		double logScore=0;
		// iterate over bases and sum the log scores.
		for (int i=0; i<bases.size(); i++)
			logScore+=Math.log10(getLikelihoodHomozygote(refAllele, altAllele, bases.get(i), qualities.get(i), genotypeQuality, maximumObservationProbability));
		return logScore;
		
		
	}

	/**
	 * The single observation version to get the homozygous likelihood.
	 * @param refAllele
	 * @param altAllele
	 * @param base
	 * @param quality
	 * @return
	 */
	private double getLikelihoodHomozygote (final byte refAllele, final byte altAllele, final Byte base, final Byte quality, final Double genotypeQuality, final Double maximumObservationProbability) {
		double prob = phredScoreToErrorProbability(quality);
		prob=getMaxErrorScore(prob, maximumObservationProbability);
		double score;
		if (base==refAllele)
			score = 1-prob;
		else
			score = prob;
		if (genotypeQuality!=null)
			score*=genotypeQuality;
		return score;
	}


	/**
	 * The one or more observation version to get the heterozygous log likelihood.
	 * Iteratively calculates the heterozygous log likelihood, takes the log, and sums those logs.
	 * @param refAllele
	 * @param altAllele
	 * @param bases
	 * @param qualities
	 * @return
	 */
	private double getLogLikelihoodHeterozygote (final byte refAllele, final byte altAllele, final List<Byte> bases, final List<Byte> qualities, final Double genotypeQuality, final Double maximumObservationProbability) {
		if (refAllele==altAllele)
			throw new IllegalArgumentException("For heterozygous likelihood, ref allele [" + refAllele +"] and alt allele [" + altAllele+ "] must not match!");

		double logScore=0;
		// iterate over bases and sum the log scores.
		for (int i=0; i<bases.size(); i++) {
			//accumulate
			double score = getLikelihoodHeterozygote(refAllele, altAllele, bases.get(i), qualities.get(i), genotypeQuality, maximumObservationProbability);
			logScore+=Math.log10(score);
		}
		return logScore;
	}

	/**
	 * The single observation version to get the heterozygous likelihood.
	 * @param refAllele
	 * @param altAllele
	 * @param base
	 * @param quality
	 * @return
	 */
	private double getLikelihoodHeterozygote (final byte refAllele, final byte altAllele, final Byte base, final Byte quality, final Double genotypeQuality, final Double maximumObservationProbability) {
		if (refAllele==altAllele)
			throw new IllegalArgumentException("For heterozygous likelihood, ref allele [" + refAllele +"] and alt allele [" + altAllele+ "] must not match!");

		double errorProb = phredScoreToErrorProbability(quality);
		errorProb=getMaxErrorScore(errorProb, maximumObservationProbability);
		double score;
		if (base==refAllele || base==altAllele)
			score=((1-errorProb)/2)+(errorProb/2);
		else
			score=errorProb;
		if (genotypeQuality!=null)
			score*=genotypeQuality;
		return (score);
	}

	private double getMaxErrorScore(final double errorProb, final Double maximumObservationProbability) {
		if (maximumObservationProbability==null)
			return errorProb;
		return Math.max(errorProb, maximumObservationProbability);
	}

	//convert to regular probability from Phred base quality: 10^(quality/-10)
	/**
	 * Converts from the phread score to the probability the base was called incorrectly.
	 * @param phreadScore
	 * @return
	 */
	public double phredScoreToErrorProbability (final byte phreadScore) {
		if (phreadScore<0 || phreadScore>127)
			throw new IllegalArgumentException("Phred score must be between 0 and 127.");
		// check cache.
		double cachedResult = this.phredToProbability[phreadScore];
		if (cachedResult!=-1)
			return cachedResult;

		// if you divide by -10 as an integer you get a rounded number.
		double d1= phreadScore/ (double)-10;
		double r = Math.pow(10, d1);
		this.phredToProbability[phreadScore]=r;
		return (r);
	}




	/**
	 * Convert from the probability that the base is an error to Phread score.
	 * a base that is P=0.9 correct (0.1 error) has a score of 10.
	 * -10*(log10(0.1))=10
	 * This value is rounded to the nearest whole number.
	 * @param errorProbability that the base is an error.
	 * @return
	 */
	public byte errorProbabilityToPhredScore (final double errorProbability) {
		if (errorProbability<0 || errorProbability>1)
			throw new IllegalArgumentException("Probability must be between 0 and 1, inclusive" );
		return (byte) Math.round((-10*(Math.log10(errorProbability))));
	}
}
