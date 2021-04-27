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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.stat.StatUtils;

import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;

/**
 * Helper methods for calculating single donor, doublet, and missing data likelihoods and converting to/from Phred scores.
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
	 * @param allLikelihoods A list of likelihoods
	 * @return the probability of the observation as defined by the max(likelihood)/sum(likelihood)
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
	 * Calculate the error rates to observe the reference and alternate allele for a single UMI/SNP/Cell
	 * 
	 * An allele may be observed from three sources:
	 * 1) The donor emits the allele
	 * 2) The sequencing chemistry makes a mistake
	 * 3) The allele comes from a contaminating transcript that was incorporated into the cell from cell free RNA
	 * 
	 * The sequencing chemistry error rate is defined by the base quality.
	 * The likelihood the read came from contamination is a combination of two factors: how common the allele is in the population, and
	 * what fraction of transcripts come from cell free RNA.
	 * 
	 * F = frequency of the minor allele in the population
	 * C = fraction of transcripts in this cell from cell free RNA
	 * B = base error rate
	 * 
	 * The alternate allele error rate is: ((F*C)+B) - (F*C*B)
	 * The reference base error rate uses the frequency of the reference allele [1-minor allele frequency].
	 * The reference allele error rate is: (((1-F)*C)+B) - (((1-F))*B*C) 
	 * 
	 * There's a numerical issue where F=1 and C=1.  In this case the error rate will also be 1.
	 * This is problematic because the likelihood is 1-error rate, and we shouldn't have 0 probabilities (especially as they are logged.)
	 * 
	 * Cap the maximum error rate at 1-(base error rate).
	 * 
	 * @param baseQuality The base quality at this position
	 * @param minorAlleleFrequency The minor allele frequency of this variant in the population
	 * @param contamination The fraction of transcripts estimated to come from cell free RNA
	 * @return A double array containing the reference and alternate allele error rates
	 */
	public double [] getContaminationErrorRates (byte baseQuality, double minorAlleleFrequency, double contamination) {
		
		double baseErrorRate = phredScoreToErrorProbability(baseQuality);		
		double maxErrorRate=1-baseErrorRate;		
		double refAlleleFrequency=1-minorAlleleFrequency;
		
		double refErrorRate=((refAlleleFrequency*contamination)+baseErrorRate) - (refAlleleFrequency*contamination*baseErrorRate);
		if (refErrorRate>(maxErrorRate)) refErrorRate=maxErrorRate;
		
		double altErrorRate=((minorAlleleFrequency*contamination)+baseErrorRate) - (minorAlleleFrequency*contamination*baseErrorRate);
		if (altErrorRate>maxErrorRate) altErrorRate=maxErrorRate;
		
		double [] result = {refErrorRate,altErrorRate};
		return (result);
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
	 * @param quality The phred scores (qualities) of the observed bases
	 * @param qualities The qualities of the bases observed.  This list must be in the same order as the bases observed.
	 * @param missingDataPenality a pre-computed likelihood to use instead of computing a likelihood for a genotype.
	 * @param genotypeProbability The probability of the genotype being correct (as set by the GQ field of the genotype in the VCF).  Can be set null to ignore.
	 * @param maximumObservationProbability If set, this is the maximum penalty that can be generated for a single observation.  Optional (set to null to ignore) 
	 * @return The likelihood of observing this set of UMIs. 
	 */
	public double getLogLikelihoodMixedModel (final char refAllele, final char altAllele, final List<GenotypeType> genotypes,
			final List<Double> mixture, final List<Byte> bases, final List<Byte> qualities, final Double missingDataPenality, final Double genotypeProbability, final Double maximumObservationProbability) {

		byte ref= StringUtil.charToByte(refAllele);
		byte alt= StringUtil.charToByte(altAllele);
		double result = getLogLikelihoodMixedModel(ref, alt, genotypes, mixture, bases, qualities, missingDataPenality, genotypeProbability, maximumObservationProbability);
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
	 * @param qualities The phred qualities of the bases observed.  This list must be in the same order as the bases observed.
	 * @param missingDataPenality a pre-computed likelihood to use instead of computing a likelihood for a genotype.
	 * @param genotypeProbability The probability of the genotype being correct (as set by the GQ field of the genotype in the VCF).  Can be set null to ignore.
	 * @param maximumObservationProbability If set, this is the maximum penalty that can be generated for a single observation.  Optional (set to null to ignore) 
	 * @return The likelihood of observing this set of UMIs. 
	 */
	public double getLogLikelihoodMixedModel (final byte refAllele, final byte altAllele, final List<GenotypeType> genotypes,
			final List<Double> mixture, final List<Byte> bases, final List<Byte> qualities, final Double missingDataPenality, final Double genotypeProbability, final Double maximumObservationProbability) {

		if (genotypes.size()!=mixture.size())
			throw new IllegalArgumentException("Genotype list and mixture list must be the same size.");

		double result=0;

		for (int i=0; i<bases.size(); i++) {
			//TODO add contamination information per cell and minor allele frequency data to integrate into the model.
			double likelihood = getLikelihoodMixedModel(refAllele, altAllele, genotypes, mixture, bases.get(i), qualities.get(i), missingDataPenality, genotypeProbability, maximumObservationProbability, null, null);			
			result+=Math.log10(likelihood);						
		}
		return result;
	}
	
	/**
	 * For a pair of donor likelihoods across many UMIs at a given mixture, compute the final likelihood of the joint model.
	 * 
	 * This is useful for optimization of many likelihoods to find the optimal mixture coefficent that produces the maximum likelihood of the joint observations.
	 *  
	 * @param likelihoods precomputed likelihoods for each donor across all UMIs.  Each row contains the likelihood of a single UMI.  Each column contains results for one donor
	 * @param mixture The mixture of the donors.  The list is the same length as the number of columns in the likelihoods argument.  This is usually expressed as
	 * a vector that sums to 1.
	 * @return The likelihood across all UMIs for donors at these mixing proportions.
	 */
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
	 * Calculates the likelihood for a list of genotype states for a single UMI observation.
	 * Here, the likelihoods for the genotype states have been precomputed.  The number of likelihoods should be equal to the number of mixture coefficients.
	 * This is the equivalent of the weighted average of the likelihoods.
	 * @param likelihoods An array of doubles representing the likelihoods of the genotypes in this mixture
	 * @param mixture A list of mixture coefficients to weight the genotype likelihoods.
	 * @return The likelihood of the mixed model.
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
	 * Calculates the likelihood for a list of genotype states for a single observation.
	 * @param ref the reference allele for the variant
	 * @param alt the alternate allele for the variant
	 * @param genotypes The list of genotypes to calculate likelihoods for
	 * @param mixture A list of mixture coefficients for the list of genotypes.
	 * @param base The observed base in the sequence pileup
	 * @param quality The phred quality of the observed base
	 * @param missingDataPenality a pre-computed likelihood to use instead of computing a likelihood for a genotype.
	 * @param genotypeProbability The likelihood of the genotype.  Can be null to ignore.
	 * @param maximumObservationProbability If set, this is the maximum penalty that can be generated for a single observation.
 	 * @param minorAlleleFrequency The minor allele frequency in the population. Optional. Set null to ignore.
	 * @param contamination The contamination rate in the population . Optional. Set null to ignore.
	 * @return the likelihood of the base/quality, given a mixture of genotypes.
	 */
	private double getLikelihoodMixedModel(final byte ref, final byte alt, final List<GenotypeType> genotypes, final List<Double> mixture,
			final byte base, final byte quality, final Double missingDataPenality, final Double genotypeProbability, final Double maximumObservationProbability, 
			Double minorAlleleFrequency, Double contamination) {

		double result=0;
		double sumMixture=0;
		// for each genotype model, calculate the likelihood, multiply by the mixture for that model, then add to the result.
		for (int i=0; i<genotypes.size(); i++) {
			GenotypeType genotype = genotypes.get(i);
			if ((genotype==GenotypeType.NO_CALL || genotype==GenotypeType.UNAVAILABLE) && missingDataPenality==null)
				throw new IllegalArgumentException("If using NO_CALL or UNAVAILABLE genotypes, must set a missingDataPenality!");
			double mix = mixture.get(i);
			sumMixture+=mix;			
			double likelihood=getLikelihood (ref, alt, genotype, base, quality, missingDataPenality, genotypeProbability, maximumObservationProbability, minorAlleleFrequency, contamination);
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
	 * Get a list of likelihoods for a given base/quality for a list of genotypes.
	 * 
	 * If the referenceAllele, minorAlleleFrequency and contamination parameters are all set, take contamination into account when calculating the error rate to see alleles.
	 * 
	 * @param ref the reference allele for the variant
	 * @param alt the alternate allele for the variant
	 * @param genotypes The list of genotypes to calculate likelihoods for
	 * @param base The observed base in the sequence pileup
	 * @param quality The phred quality of the observed base
	 * @param missingDataPenality a pre-computed likelihood to use instead of computing a likelihood for a genotype.
	 * @param genotypeProbability The likelihood of the genotype.  Can be null to ignore.
	 * @param maximumObservationProbability If set, this is the maximum penalty that can be generated for a single observation.
	 * @param minorAlleleFrequency The minor allele frequency in the population. Optional. Set null to ignore.
	 * @param contamination The contamination rate in the population . Optional. Set null to ignore.
	 * @return An array of likelihoods in the same order as the submitted genotypes.
     */
	public double [] getLikelihoodManyObservations (final byte ref, final byte alt, final List<GenotypeType> genotypes, final byte base, final byte quality, 
			final Double missingDataPenality, final Double genotypeProbability, final Double maximumObservationProbability,
			Double minorAlleleFrequency, Double contamination) {
		
		double [] result = new double [genotypes.size()];
		// this doesn't implement caching as the typical use case is a genotype list of 2 different genotypes for doublet detection.
		for (int i=0; i<genotypes.size(); i++) {
			GenotypeType genotype = genotypes.get(i);
			result[i]=getLikelihood (ref, alt, genotype, base, quality, missingDataPenality, genotypeProbability, maximumObservationProbability,
					minorAlleleFrequency, contamination);			
		}				
		return result;		
	}
	
	

	/**
	 * Calculates the likelihood for a single donor genotype for one or more UMI observations.
	 * 
	 * If the referenceAllele, minorAlleleFrequency and contamination parameters are all set, take contamination into account when calculating the error rate to see alleles.
	 * @param alleleOne The first allele of the SNP for the sample
	 * @param alleleTwo The second allele of the SNP for the sample
	 * @param bases The observed bases in the sequence pileup
	 * @param qualities	The phred quality of the observed base
	 * @param genotypeProbability The likelihood of the genotype.  Can be null to ignore.
	 * @param maximumObservationProbability If set, this is the maximum penalty that can be generated for a single observation.  This is ignored if the contamination parameters are set.
	 * @param referenceAllele The reference allele for this variant in the population. Optional. Set null to ignore.
	 * @param minorAlleleFrequency The minor allele frequency in the population. Optional. Set null to ignore.
	 * @param contamination The contamination rate in the population . Optional. Set null to ignore.
  
	 * @return the likelihood of the bases/qualities, given the alleles observed.
	 */
	public double getLogLikelihood (final char alleleOne, final char alleleTwo, final List<Byte> bases, final List<Byte> qualities, final Double genotypeProbability, final Double maximumObservationProbability, 
			final Byte referenceAllele, Double minorAlleleFrequency, Double contamination) {
		byte ref= StringUtil.charToByte(alleleOne);
		byte alt= StringUtil.charToByte(alleleTwo);
		return getLogLikelihood(ref, alt, bases, qualities, genotypeProbability, maximumObservationProbability, referenceAllele, minorAlleleFrequency, contamination);
	}

	/**
	 * Calculates the likelihood for a single donor genotype for one or more UMI observations.
	 * 
	 * This dispatches to the proper likelihood calculation depending on which parameters are not null.
	 * 
	 * If the referenceAllele, minorAlleleFrequency and contamination parameters are all set, take contamination into account when calculating the error rate to see alleles.
	 * 
	 * 
	 * @param alleleOne The first allele of the SNP for the sample
	 * @param alleleTwo The second allele of the SNP for the sample
	 * @param bases The observed bases in the sequence pileup
	 * @param qualities	The phred quality of the observed base
	 * @param genotypeProbability The likelihood of the genotype.  Optional. Set null to ignore.
	 * @param maximumObservationProbability If set, this is the maximum penalty that can be generated for a single observation.  This is ignored if the contamination parameters are set.
	 * @param referenceAllele The reference allele for this variant in the population. Optional. Set null to ignore.
	 * @param minorAlleleFrequency The minor allele frequency in the population. Optional. Set null to ignore.
	 * @param contamination The contamination rate in the population . Optional. Set null to ignore.
	 * @return the likelihood of the bases/qualities, given the alleles observed.
	 */
	public double getLogLikelihood (final byte alleleOne, final byte alleleTwo, final List<Byte> bases, final List<Byte> qualities, final Double genotypeProbability, final Double maximumObservationProbability, 
			final Byte referenceAllele, Double minorAlleleFrequency, Double contamination) {
		double logScore=0;
		// iterate over bases and sum the log scores.
		for (int i=0; i<bases.size(); i++)
			logScore+=Math.log10(getLikelihood(alleleOne, alleleTwo, bases.get(i), qualities.get(i), genotypeProbability, maximumObservationProbability, referenceAllele, minorAlleleFrequency, contamination));
		return logScore;							
	}
	
	
	/**
	 * Convenience method to get a likelihood for a given base/quality for a genotype.
	 * 
	 * This is parameterized to run on a Genotype object in combination with the reference and alternate allele.
	 * 
	 * If the referenceAllele, minorAlleleFrequency and contamination parameters are all set, take contamination into account when calculating the error rate to see alleles.
	 * 
	 * @param ref the reference allele for the variant
	 * @param alt the alternate allele for the variant
	 * @param genotype The genotype to calculate likelihoods for
	 * @param base The observed base in the sequence pileup
	 * @param quality The phred quality of the observed base
	 * @param missingDataPenality a pre-computed likelihood to use instead of computing a likelihood for a genotype.
	 * @param genotypeProbability The likelihood of the genotype.  Can be null to ignore.
	 * @param maximumObservationProbability If set, this is the maximum penalty that can be generated for a single observation.
	 * @param minorAlleleFrequency The minor allele frequency in the population. Optional. Set null to ignore.
	 * @param contamination The contamination rate in the population . Optional. Set null to ignore.
	 * @return An array of likelihoods in the same order as the submitted genotypes.
     */
	public double getLikelihood (final byte ref, final byte alt, final GenotypeType genotype, final byte base, final byte quality, 
			final Double missingDataPenality, final Double genotypeProbability, final Double maximumObservationProbability,
			Double minorAlleleFrequency, Double contamination) {
		
		if ((genotype==GenotypeType.NO_CALL || genotype==GenotypeType.UNAVAILABLE) && missingDataPenality==null)
			throw new IllegalArgumentException("If using NO_CALL or UNAVAILABLE genotypes, must set a missingDataPenality!");
		double likelihood = 0;
		// calculate the single likelihood.
		switch (genotype) {
		case HOM_REF:
			likelihood=getLikelihood(ref, ref, base, quality, genotypeProbability, maximumObservationProbability, ref, minorAlleleFrequency, contamination); break;				
		case HOM_VAR:
			likelihood=getLikelihood(alt, alt, base, quality, genotypeProbability, maximumObservationProbability, ref, minorAlleleFrequency, contamination); break;				
		case HET:
			likelihood=getLikelihood(ref, alt, base, quality, genotypeProbability, maximumObservationProbability, ref, minorAlleleFrequency, contamination); break;
		case NO_CALL:
			likelihood=missingDataPenality; break;
		case UNAVAILABLE:
			likelihood=missingDataPenality; break;
		default:
		}
		return (likelihood);
	}
	
	/**
	 * Convenience method to dispatch likelihood calculation for a single UMI to the correct likelihood method. 
	 * 
	 * @param alleleOne
	 * @param alleleTwo
	 * @param base
	 * @param quality
	 * @param genotypeProbability
	 * @param maximumObservationProbability
	 * @param referenceAllele
	 * @param minorAlleleFrequency
	 * @param contamination
	 * @return
	 */
	public double getLikelihood (final byte alleleOne, final byte alleleTwo, final byte base, final byte quality, final Double genotypeProbability, final Double maximumObservationProbability, 
			final Byte referenceAllele, Double minorAlleleFrequency, Double contamination) {
		if (alleleOne==alleleTwo) {
			if (contamination!=null && minorAlleleFrequency!=null && referenceAllele!=null) 
				return getLikelihoodHomozygoteWithContamination(alleleOne, alleleTwo, base, quality, genotypeProbability, referenceAllele, minorAlleleFrequency, contamination);			
			return getLikelihoodHomozygote(alleleOne, alleleTwo, base, quality, genotypeProbability, maximumObservationProbability);	
		}			
		else {
			if (contamination!=null && minorAlleleFrequency!=null && referenceAllele!=null)
				return getLikelihoodHeterozygoteWithContamination(alleleOne, alleleTwo, base, quality, genotypeProbability, referenceAllele, minorAlleleFrequency, contamination);
			return getLikelihoodHeterozygote(alleleOne, alleleTwo, base, quality, genotypeProbability, maximumObservationProbability);
		}		
	}

	
	
	/****************************************************************************************************
	 * Likelihood calculations for homozygous and heterozygous sites
	 * Single observations -> likelihood or many observations summarized by the sum log10(likelihoods)
	 * With or without factoring in contamination
	 ****************************************************************************************************/
	
	/**
	 * Calculate the likelihood of a single UMI at a homozygous genotype
	 * 
	 * By default, this uses the base error to determine the likelihood - if a base has an error rate of 0.01, then when the base and genotype agree the likelihood is 1- error rate,
	 * otherwise the likelihood = the error rate.
	 *  
	 * @param alleleOne The first allele of the SNP for the sample
	 * @param alleleTwo The second allele of the SNP for the sample
	 * @param base The observed base
	 * @param quality A list of observed phred quality for the base.
	 * @param genotypeProbability The confidence in the genotype.  If set to null this is ignored.
	 * @param maximumObservationProbability A cap on the maximum likelihood penalty at any one UMI. If set to null this is ignored.  
	 * @return The likelihood for one observation
	 */
	private double getLikelihoodHomozygote (final byte alleleOne, final byte alleleTwo, final Byte base, final Byte quality, final Double genotypeProbability, final Double maximumObservationProbability) {
		double errrorProb = phredScoreToErrorProbability(quality);
		errrorProb=getMaxErrorScore(errrorProb, maximumObservationProbability);
		return (getLikelihoodHomozygote(alleleOne, alleleTwo, base, errrorProb, genotypeProbability));
	}
	
	/**
	 * Compute the likelihood of a homozygous allele, and take into account contamination with cell free RNA.
	 * This likelihood depends on both the error rate of the base AND the contamination of this cell and allele frequency data of this SNP.
	 * 
	 * The penalty score will change for every UMI (different base qualities), 
	 * every SNP (different allele frequency), 
	 * and every cell (different contamination rates)
	 * 
	 * @param alleleOne The first allele of the SNP for the sample
	 * @param alleleTwo The second allele of the SNP for the sample
	 * @param base The observed base
	 * @param quality A list of observed phred quality for the base.
	 * @param genotypeProbability The confidence in the genotype.  If set to null this is ignored.
	 * @param referenceAllele The reference allele for this variant in the population.
	 * @param minorAlleleFrequency The minor allele frequency in the population
	 * @param contamination The contamination rate in the population 
	 * @return @see getContaminationErrorRates
	 */
	double getLikelihoodHomozygoteWithContamination (final byte alleleOne, final byte alleleTwo, final byte base, final byte quality, final Double genotypeProbability, 
			final byte referenceAllele, double minorAlleleFrequency, double contamination) {
		
		double [] penaltyScores = getContaminationErrorRates(quality, minorAlleleFrequency, contamination);
		// If this donor has the reference allele, use the referencee allele penalty.
		double errrorProb = penaltyScores[0];
		// otherwise use the alternate allele penalty
		if (base!=referenceAllele) errrorProb = penaltyScores[1];		
		double result = getLikelihoodHomozygote(alleleOne, alleleTwo, base, errrorProb, genotypeProbability);
		return (result);
	}
	
	/**
	 * Simple computation of the homozygous likelihood of a single observation.
	 * If the allele matches the sequencing base return 1-error probability, else return the error probability.
	 * 
	 * @param alleleOne The first allele of the SNP for the sample
	 * @param alleleTwo The second allele of the SNP for the sample
	 * @param base The observed base in the sequencing data
	 * @param errorPenalty The likelihood score when the donor allele doesn't match the observed allele in the sequencing data.
	 * @param genotypeProbability The confidence in the genotype.  If set to null this is ignored.
	 * @return
	 */
	private double getLikelihoodHomozygote (final byte alleleOne, final byte alleleTwo, final byte base, final double errorPenalty, final Double genotypeProbability) {
		if (alleleOne!=alleleTwo)
			throw new IllegalArgumentException("For homozygous likelihood, ref allele [" + alleleOne +"] and alt allele [" + alleleTwo+ "] must match!");
		double score;
		if (base==alleleOne)
			score = 1-errorPenalty;
		else
			score = errorPenalty;
		if (genotypeProbability!=null)
			score*=genotypeProbability;
		return score;
	}
		
	/**
	 * Calculate the likelihood of a single UMI at a heterozygous genotype
	 * 
	 * By default, this uses the base error to determine the likelihood - if a base has an error rate of 0.01, then when the base and genotype agree the likelihood is 1- error rate,
	 * otherwise the likelihood = the error rate.
	 * 
	 * @param alleleOne The first allele of the SNP for the sample
	 * @param alleleTwo The second allele of the SNP for the sample
	 * @param base The observed base
	 * @param quality A list of observed phred quality for the base.
	 * @param genotypeProbability The confidence in the genotype.  If set to null this is ignored.
	 * @param maximumObservationProbability A cap on the maximum likelihood penalty at any one UMI. If set to null this is ignored.  
	 * @return The likelihood for one observation
	 */
	private double getLikelihoodHeterozygote (final byte alleleOne, final byte alleleTwo, final Byte base, final Byte quality, final Double genotypeProbability, final Double maximumObservationProbability) {
		if (alleleOne==alleleTwo)
			throw new IllegalArgumentException("For heterozygous likelihood, ref allele [" + alleleOne +"] and alt allele [" + alleleTwo+ "] must not match!");

		double errorProb = phredScoreToErrorProbability(quality);
		errorProb=getMaxErrorScore(errorProb, maximumObservationProbability);
		double score;
		if (base==alleleOne || base==alleleTwo)
			score=((1-errorProb)/2)+(errorProb/2);
		else
			score=errorProb;
		if (genotypeProbability!=null)
			score*=genotypeProbability;
		return (score);
	}
	
	/**
	 * Compute the likelihood of a heterozygous allele, and take into account contamination with cell free RNA.
	 * This likelihood depends on both the error rate of the base AND the contamination of this cell and allele frequency data of this SNP.
	 * 
	 * @param alleleOne The first allele of the SNP for the sample
	 * @param alleleTwo The second allele of the SNP for the sample
	 * @param base The observed base
	 * @param quality A list of observed phred quality for the base.
	 * @param genotypeProbability The confidence in the genotype.  If set to null this is ignored.
	 * @param maximumObservationProbability A cap on the maximum likelihood penalty at any one UMI. If set to null this is ignored.  
	 * @return The likelihood for one observation
	 */
	double getLikelihoodHeterozygoteWithContamination (final byte alleleOne, final byte alleleTwo, final byte base, final byte quality, final Double genotypeProbability, 
			final byte referenceAllele, double minorAlleleFrequency, double contamination) {
		
		if (alleleOne==alleleTwo)
			throw new IllegalArgumentException("For heterozygous likelihood, ref allele [" + alleleOne +"] and alt allele [" + alleleTwo+ "] must not match!");

		double ref;
		double alt;
		if (base==alleleOne) {
			ref=getLikelihoodHomozygoteWithContamination(alleleOne, alleleOne, alleleOne, quality, genotypeProbability, referenceAllele, minorAlleleFrequency, contamination);
			alt=getLikelihoodHomozygoteWithContamination(alleleOne, alleleOne, alleleTwo, quality, genotypeProbability, referenceAllele, minorAlleleFrequency, contamination);				
		} else {
			ref=getLikelihoodHomozygoteWithContamination(alleleTwo, alleleTwo, alleleOne, quality, genotypeProbability, referenceAllele, minorAlleleFrequency, contamination);
			alt=getLikelihoodHomozygoteWithContamination(alleleTwo, alleleTwo, alleleTwo, quality, genotypeProbability, referenceAllele, minorAlleleFrequency, contamination);
		}
		
		// TODO: UNIT TEST THIS: 
		// When at a heterozygous site, the likelihood of the donor should be higher when observing the allele that is less frequent in the population.
		// That likelihood difference should be encoded in the error rates of each allele.
		double score;
		if (base==alleleOne || base==alleleTwo) 			
			score=(ref+alt)/2;															
		else  {
			// otherwise pick the bigger penalty, giving the lowest likelihood. 
			score= Math.min(ref, alt);			
		}			
		if (genotypeProbability!=null)
			score*=genotypeProbability;
		return (score);
	}

	/**
	 * Convenience method to cap the maximum error penalty for a single UMI at a threshold.
	 * @param errorProb The original error probability
	 * @param maximumObservationProbability The maximum error probability
	 * @return The maximum value of the error probability and the given maximum
	 */
	private double getMaxErrorScore(final double errorProb, final Double maximumObservationProbability) {
		if (maximumObservationProbability==null)
			return errorProb;
		return Math.max(errorProb, maximumObservationProbability);
	}

	//convert to regular probability from Phred base quality: 10^(quality/-10)
	/**
	 * Converts from the phread score to the probability the base was called incorrectly.
	 * @param phreadScore The phred score
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
