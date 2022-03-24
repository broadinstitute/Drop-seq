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

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.GenotypeType;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.LikelihoodUtils;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPIntervalRecordI;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPUMIBasePileup;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SummarizeUMIBaseQualities;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.editdistance.MapBarcodesByEditDistance;
import picard.annotation.LocusFunction;

import java.util.*;

/**
 * Holds the bases and base qualities at a SNP position for a cell.
 * This holds the collection of pileup information across all molecular barcodes.
 *
 * @author nemesh
 *
 */
public class SampleGenotypeProbabilities implements SNPIntervalRecordI {

	private final Interval snpInterval;
	private final String cell;
	private List<SNPUMIBasePileup> pileups;
	private List<Byte> bases;
	private List<Byte> qualities;
	private boolean isModified;
	// private LikelihoodCache cachedLikelihood;
 
	public SampleGenotypeProbabilities(final Interval snpInterval, final String cell) {
		this.snpInterval=snpInterval;
		this.cell=cell;
		this.pileups=new ArrayList<>();
		// this.cachedLikelihood=new LikelihoodCache();
	}

	public List<SNPUMIBasePileup> getBackingPileups () {
		return this.pileups;
	}

	public void add (final SNPUMIBasePileup p) {
		int numBases = p.getNumBases();
		if (numBases==0)
			return; // don't add empty pileups.
		if (!p.getSNPInterval().equals(this.snpInterval)  || !p.getCell().equals(this.cell))
			throw new IllegalArgumentException("Can't add pileup " + p +" to this probabilities collection " + this);
		isModified=true;
		pileups.add(p);
	}

	public Interval getSNPInterval () {
		return this.snpInterval;
	}

	public String getCell() {
		return cell;
	}

	public List<Byte> getBases() {
		if (isModified) recalculateBasePileup();
		return bases;
	}
	
	@Override
	public int getNumBases() {
		if (isModified) recalculateBasePileup();
		return bases.size();
	}

	/**
	 * Get the number of umi observations on each base.
	 */
	public ObjectCounter<Character> getUMIBaseCounts() {
		List<Byte> bases = getBases();
		ObjectCounter<Character> result = new ObjectCounter<>();
		for (Byte b: bases)
			result.increment(StringUtil.byteToChar(b));
		return result;
	}

	/**
	 * Get the number of read observations on each base.
	 *
	 */
	public ObjectCounter<Character> getReadBaseCounts() {
		ObjectCounter<Character> result = new ObjectCounter<>();
		for (SNPUMIBasePileup p:  this.pileups) {
			List<Character> bases = p.getBasesAsCharacters();
			for (Character c: bases)
				result.increment(c);
		}
		return result;
	}

	public List<Byte> getQualities() {
		if (isModified) recalculateBasePileup();
		return this.qualities;
	}

	/**
	 * From the given set of pileup objects, construct the final list of consensus bases and qualities for those pileups.
	 * There will be one base/quality score per pileup.
	 * This should be triggered whenever functions call for bases or quality scores.
	 */
	private void recalculateBasePileup() {
		this.bases=new ArrayList<>();
		this.qualities=new ArrayList<>();

		for (SNPUMIBasePileup p: this.pileups) {
			SummarizeUMIBaseQualities summarize = new SummarizeUMIBaseQualities(p.getBases(), p.getQualities());
			byte base = summarize.getMostCommonBase();
			byte qual = (byte) summarize.getSummarizedPhredScore();
			bases.add(base);
			qualities.add(qual);
		}

		this.isModified=false;
	}
		
	/**
	 * For this pileup, emit the likelihood that the reads were drawn from a diploid genotype.
	 * Used the fixed error rate if it is non-null, otherwise use the base qualities for error rates.
	 * Multiply the likelihood by the genotypeQuality, once it has been converted from a phred based score to a probability between 0 and 1.
	 * 
	 * If you wish to model cell free RNAs effect on likelihood, the referenceAllele, minorAlleleFrequency, and contamination must all be not null.
	 * If these parameters are set, then the maximumObservationProbability parameter is ignored even if set.
	 * 
	 * @param alleleOne The first allele of the SNP for the sample
	 * @param alleleTwo The second allele of the SNP for the sample
	 * @param fixedGenotypeErrorRate Instead of using the base qualities to compute likelihoods, used a fixed error rate instead.  Optional (set to null to ignore)
	 * @param maximumObservationProbability If set, this is the maximum penalty that can be generated for a single observation.  Optional (set to null to ignore)
	 * @param referenceAllele The reference allele for this variant in the population. (Optional, set null to ignore).   
	 * @param minorAlleleFrequency The minor allele frequency in the population (Optional, set null to ignore).
	 * @param contamination The fraction of UMIs estimated to come from ambient RNA for this cell. (Optional, set null to ignore).

	 * @return The likelihood of the genotype given the data.
	 */
	public double getLogLikelihood (final char alleleOne, final char alleleTwo, final Double fixedGenotypeErrorRate, final Double genotypeProbability, 
			final Double maximumObservationProbability, final Character referenceAllele, Double minorAlleleFrequency, Double contamination) {
		List<Byte> quals;
		if (fixedGenotypeErrorRate!=null) {
			byte phreadScore = LikelihoodUtils.getInstance().errorProbabilityToPhredScore(fixedGenotypeErrorRate);
			Byte [] tempQuals = new Byte [getBases().size()];
			Arrays.fill(tempQuals, phreadScore);
			quals=Arrays.asList(tempQuals);
		} else
			quals=getQualities();
		System.out.println(String.format("%c %c %f %f %f %c %f %f {%s} {%s}", alleleOne, alleleTwo, fixedGenotypeErrorRate,
				genotypeProbability, maximumObservationProbability, referenceAllele, minorAlleleFrequency, contamination,
				StringUtil.join(",", getBases()), StringUtil.join(",", quals)));
		//TODO: add in contamination parameters
		double likelihood = LikelihoodUtils.getInstance().getLogLikelihood(alleleOne, alleleTwo, getBases(), quals,
				genotypeProbability, maximumObservationProbability, referenceAllele, minorAlleleFrequency, contamination);
		return likelihood;
	}
			
	/**
	 * For this pileup, emit the missing data likelihood
	 * This uses a list of genotype states and a list of mixture ratios for those states.
	 * @param alleleOne The first allele of the SNP for the sample
	 * @param alleleTwo The second allele of the SNP for the sample
	 * @param genotypes The list of genotype states for the observed genotypes for the variant (donors that correctly genotyped)
	 * @param mixture The list of mixture parameters that these genotypes are mixed by.  This list must be in the same order as the genotype states.
	 * @param genotypeProbability the probability of the genotype between 0 and 1, inclusive.  Can be set to null to ignore.
	 * @param maximumObservationProbability If set, this is the maximum penalty that can be generated for a single observation.  Optional (set to null to ignore)
	 * @param minorAlleleFrequency The minor allele frequency in the population. Optional. Set null to ignore.
	 * @param contamination The contamination rate in the population . Optional. Set null to ignore.
	 * @return the likelihood for any donors with missing data.
	 */
	public double getLogLikelihoodMissingData (final char alleleOne, final char alleleTwo, final List<GenotypeType> genotypes,
			final List<Double> mixture, final Double fixedGenotypeErrorRate, final Double genotypeProbability,  
			final Double maximumObservationProbability, Double minorAlleleFrequency, Double contamination) {

		List<Byte> quals;
		if (fixedGenotypeErrorRate!=null) {
			byte phreadScore = LikelihoodUtils.getInstance().errorProbabilityToPhredScore(fixedGenotypeErrorRate);
			Byte [] tempQuals = new Byte [getBases().size()];
			Arrays.fill(tempQuals, phreadScore);
			quals=Arrays.asList(tempQuals);
		} else
			quals=getQualities();

		return LikelihoodUtils.getInstance().getLogLikelihoodMixedModel(alleleOne, alleleTwo, genotypes, mixture, getBases(), quals, fixedGenotypeErrorRate, 
				genotypeProbability, maximumObservationProbability, minorAlleleFrequency, contamination );
	}



	public List<Character> getBasesAsCharacters() {
		if (isModified) recalculateBasePileup();
		if (this.bases==null || this.bases.size()==0) return Collections.EMPTY_LIST;

		List<Character> characters = new ArrayList<>(this.bases.size());
		for (Byte b: this.bases)
			characters.add(StringUtil.byteToChar(b));
		return characters;
	}


	/**
	 * This object consists of a collection of pileups with different molecular barcodes.
	 * Since those barcodes can be within <edit distance> of each other, it can be desirable to
	 * collapse neighboring pileups together into single data sets.
	 * @param editDistance The hamming distance to collapse pileups.
	 */
	public void collapseUMIs (final int editDistance) {
		// short circuit for ED=0.
		if (editDistance==0)
			return;

		// maps the UMI to the pileup.
		Map<String, SNPUMIBasePileup> pileMap = new HashMap<>();
		// the pileup #bases keyed by molecular barcode.
		ObjectCounter<String> umis = new ObjectCounter<>();
		for (SNPUMIBasePileup p: this.pileups) {
			pileMap.put(p.getMolecularBarcode(), p);
			umis.increment(p.getMolecularBarcode());
		}
		// do any of the UMIs collapse?  Map them here.
		MapBarcodesByEditDistance med = new MapBarcodesByEditDistance(false, 1, 0);
		Map<String, List<String>> mapping = med.collapseBarcodes(umis, false, editDistance);

		// merge UMI results - the number of final pileups is the number of mapping keys.
		List<SNPUMIBasePileup> resultPileUps = new ArrayList<>(mapping.keySet().size());

		// loop over each key, and combine the key + any values.
		for (String primaryMolBC: mapping.keySet()) {
			SNPUMIBasePileup primaryPile = pileMap.get(primaryMolBC);
			List<Byte> bases = primaryPile.getBases();
			List<Byte> quals = primaryPile.getQualities();

			// if there are additional pileups, add them to the primary
			for (String secondaryUMI: mapping.get(primaryMolBC)) {
				SNPUMIBasePileup secondaryPile = pileMap.get(secondaryUMI);
				bases.addAll(secondaryPile.getBases());
				quals.addAll(secondaryPile.getQualities());
			}
			resultPileUps.add(primaryPile);
		}

		this.pileups=resultPileUps;
	}

	public Set<LocusFunction> getLocusFunctions() {
		Set<LocusFunction> result = new HashSet<>();
		for (SNPUMIBasePileup p: this.pileups)
			result.addAll(p.getLocusFunctions());
		return result;
	}

	@Override
	public String toString () {
		if (isModified) recalculateBasePileup();
		StringBuilder b = new StringBuilder();
		b.append("snp [" + this.snpInterval.toString() +"] ");
		b.append("cell [" + this.cell +"] ");
		b.append(this.getBasesAsCharacters().toString()+ " ");
		if (this.qualities!=null) b.append(this.qualities);
		Set<LocusFunction> locusFunctionSet = getLocusFunctions();
		if (locusFunctionSet!=null && locusFunctionSet.size()>0)
			b.append("locus function " + locusFunctionSet +"");
		return b.toString();
	}
	
	/**
	 * Instead of calculating likelihoods over and over, cache the results for the possible alleles.
	 * This takes the two alleles, concatenates them into a single string, and stores the result.
	 *
	 * Turns out, this is way more expensive than just calculating the likelihoods over and over....caching is premature optimization.
	 * @author nemesh
	 *
	 */
	/*
	class LikelihoodCache {
		// have this be a byte [] of the two alleles.
		private Map<char [], Double> cache = new HashMap<char [], Double>();
		
		private Double cachedfixedGenotypeErrorRate=null;
		private Double cachedGenotypeProbability=null;
		private Double cachedMaximumObservationProbability=null;		
		
		private char [] getKey (final char refAllele, final char altAllele) {
			char [] alleles = {refAllele, altAllele};
			return alleles;
		}
		
		/**
		 * If the fixedGenotypeErrorRate, genotypeProbability and maximumObservationProbability match a previous call for a given reference and alternate allele, then return
		 * the cached likelihood score.
		 * @param refAllele Allele 1 of the genotype of the donor 
		 * @param altAllele Allele 2 of the genotype of the donor 
		 * @param bases A list of bases observed at a sequence pileup to compute the likelihood of this ref/alt allele.
		 * @param quals A list of quality scores for the bases.  
		 * @param fixedGenotypeErrorRate Instead of using the base qualities to compute likelihoods, used a fixed error rate instead.   Set to null to ignore.
		 * @param genotypeQuality the genotype quality formatted as a phred score.  Set to null to ignore.
		 * @param maximumObservationProbability The maximum penalty score that can be assigned to any one UMI observation.  Set to null to ignroe.
		 * @param maximumObservationProbability
		 * @return The likelihood of the bases at the given quals, given the ref and alt alleles of the donor
		 */
		/*
		public double getLikelihood(final char refAllele, final char altAllele, final List<Byte> bases,  List<Byte> quals, 
				final Double fixedGenotypeErrorRate, final Double genotypeProbability, final Double maximumObservationProbability) {
			// Are the settings the same for this call as a previous call?
			
			if (useCache(fixedGenotypeErrorRate, genotypeProbability, maximumObservationProbability)) {
				// retrieve from map.
				Double result = cache.get(getKey(refAllele, altAllele));
				if (result!=null) return (result);
			}
			// no cached result found, compute it up and store.
			double likelihood = LikelihoodUtils.getInstance().getLogLikelihood(refAllele, altAllele, getBases(), quals, genotypeProbability, maximumObservationProbability);
			this.cachedfixedGenotypeErrorRate=fixedGenotypeErrorRate;
			this.cachedGenotypeProbability=genotypeProbability;
			this.cachedMaximumObservationProbability=maximumObservationProbability;
			cache.put(getKey(refAllele, altAllele), likelihood);
			return likelihood;			
		}
		
		public void setLikelihood (final char refAllele, final char altAllele, double likelihood, final Double fixedGenotypeErrorRate, final Double genotypeProbability, final Double maximumObservationProbability) {
			cache.put(getKey(refAllele, altAllele), likelihood);
		}
		
		private boolean useCache (final Double fixedGenotypeErrorRate, final Double genotypeProbability, final Double maximumObservationProbability) {
			return true;
			//return (Objects.equals(fixedGenotypeErrorRate, cachedfixedGenotypeErrorRate) && Objects.equals(genotypeProbability, cachedGenotypeProbability) 
			//		&& Objects.equals(maximumObservationProbability, cachedMaximumObservationProbability));
		}
		
		
	}
	*/
}


