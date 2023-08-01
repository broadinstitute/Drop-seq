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
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.apache.commons.math3.stat.descriptive.summary.Sum;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.LikelihoodUtils;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.CellAssignmentUtils;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.SampleGenotypeProbabilities;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.GenotypeType;

public class VariantDataFactory {

	private final List<SampleGenotypeProbabilities> allProbs;
	private final GenotypeMatrix genotypeMatrix;
	private final String cell;
	private final Double fixedErrorRate;
	private final Map<SampleGenotypeProbabilities, Double> missingDataPenalties;
	private final Double maximumObservationProbability;
	private final Double contamination;
	private final Map<Interval, Double> variantMinorAlleleFrequency;
	
	
	public VariantDataFactory (final String cell, final List<SampleGenotypeProbabilities> probs, final GenotypeMatrix genotypeMatrix) {
		this(cell, probs, genotypeMatrix, null, false, null, null, null);
	}

	public VariantDataFactory (final String cell, final List<SampleGenotypeProbabilities> probs, final GenotypeMatrix genotypeMatrix, final Double fixedErrorRate) {
		this(cell, probs, genotypeMatrix, fixedErrorRate, false, null, null, null);
	}

	/**
	 * Constructs a factory that generates collections of variant data that can be optimized for a mixture of samples.
	 *
	 * @param cell The cell barcode that all the data comes from
	 * @param probs The pileups of SNP data
	 * @param genotypeMatrix A matrix of genotype states
	 * @param fixedErrorRate A fixed error rate between 0 and 1 or null.
	 * @param maximumObservationProbability cap the error rate per UMI at this value globally
	 * @param cellContaminationMap Use the estimated ambient RNA contamination in cells to modify likelihood .  Key=cell barcode, value=maximum error rate.
	 * @param variantMinorAlleleFrequency The estimated minor allele frequency of each variant.  
	 * @param useMissingDataPenalty Should genotypes that are set to GenotypeType.NO_CALL have a global (across all samples for the missing snp) value generated and used? See @getMissingDataPenalities
	 */
	public VariantDataFactory (final String cell, final List<SampleGenotypeProbabilities> probs, final GenotypeMatrix genotypeMatrix, final Double fixedErrorRate, 
			final boolean useMissingDataPenalty, final Double maximumObservationProbability, Map<String,Double> cellContaminationMap, Map<Interval, Double> variantMinorAlleleFrequency) {
		
		this.contamination = CellAssignmentUtils.getNullableValue(cellContaminationMap, cell);
		this.maximumObservationProbability=maximumObservationProbability;
		this.variantMinorAlleleFrequency = variantMinorAlleleFrequency;
		
		this.genotypeMatrix=genotypeMatrix;
		this.cell=cell;
		this.fixedErrorRate=fixedErrorRate;
		
		// initialize storage.
		this.allProbs=new ArrayList<>(probs.size());
		for (SampleGenotypeProbabilities p: probs)
			// validate that the probabilities have the cell you say you're using.
			if (p.getCell().equals(cell))
				allProbs.add(p);
			else
				throw new IllegalArgumentException("While populating data for cell " + cell+ " saw pileup data for cell " + p.getCell());
		// if we want to use missing data penalties, set all of them.
		if (useMissingDataPenalty)
			this.missingDataPenalties=getMissingDataPenalities(probs, genotypeMatrix, fixedErrorRate, maximumObservationProbability);
		else // empty map.
			this.missingDataPenalties=Collections.emptyMap();
	}

	/**
	 * Generates a collection of variants for the two samples.
	 * Only generates a list with elements if both sampleOne and sampleTwo are in the genotypeMatrix.
	 * Otherwise produces an empty list.
	 * Only add a variant if it has at least one called REF/HET/VAR.
	 * @param sampleOne
	 * @param sampleTwo
	 * @return
	 */
	public VariantDataCollection getVariantData (final String sampleOne, final String sampleTwo) {

		// if one of the genotype states is missing, don't add.
		// only test pairs where the data was set - donors
		// to be included should have at least a missing value.
		if (!this.genotypeMatrix.containsDonor(sampleOne) || !this.genotypeMatrix.containsDonor(sampleTwo)) {
			List<VariantData> empty = Collections.emptyList();
			return new VariantDataCollection(empty, sampleOne, sampleTwo, this.cell);
		}
		
		// this only works in parallel if the allele freqs in the genotype matrix are precomputed.
		List<VariantData> d= this.allProbs.parallelStream().map(x -> getVariantData(x, sampleOne, sampleTwo)).filter(x -> x!=null).collect(Collectors.toList());
		VariantDataCollection result = new VariantDataCollection(d, sampleOne, sampleTwo, this.cell);
		return result;				
	}
	
	private VariantData getVariantData (SampleGenotypeProbabilities p, final String sampleOne, final String sampleTwo) {
		Interval i = p.getSNPInterval();
		Double maf = CellAssignmentUtils.getNullableValue(this.variantMinorAlleleFrequency, i);
		
		double [] genotypeFreqs = genotypeMatrix.getGenotypeFrequencies(i);
		double sum = new Sum().evaluate(genotypeFreqs);
		// need at least one ref/het/var call in this variant to continue.
		// in practice the variants are never 
		if (sum==0) 
			return (null);
		GenotypeType s1 = genotypeMatrix.getGenotype(i, sampleOne);
		GenotypeType s2 = genotypeMatrix.getGenotype(i, sampleTwo);

		byte refAllele = genotypeMatrix.getRefBase(i);
		byte altAllele = genotypeMatrix.getAltBase(i);

		List<Byte> basesFinal = new ArrayList<>();
		List<Byte> qualsFinal = new ArrayList<>();
		List<Byte> b = p.getBases();
		List<Byte> qual = p.getQualities();
		// filter
		for (int j=0; j<b.size(); j++) {
			byte currentB=b.get(j);
			if (currentB==refAllele || currentB==altAllele) {
				basesFinal.add(currentB);
				qualsFinal.add(qual.get(j));
			}
		}

		Double missingDataPenality = missingDataPenalties.get(p);
		// should we have the possibility of both donors having a missing value?
		// that would give any two donor pairs the same number of SNPs.
		if ((s1!=GenotypeType.NO_CALL && s2!=GenotypeType.NO_CALL) || (missingDataPenality!=null)) {
			if (this.fixedErrorRate!=null) {
				byte phreadScore = LikelihoodUtils.getInstance().errorProbabilityToPhredScore(this.fixedErrorRate);
				Byte [] q = new Byte [qualsFinal.size()];
				Arrays.fill(q, phreadScore);
				qualsFinal=Arrays.asList(q);
			}
			VariantData vd =new VariantData(i, refAllele, altAllele, s1, s2, basesFinal, qualsFinal, missingDataPenality, maximumObservationProbability, this.contamination, maf);
			// only add a variant if it has data on the reference or alternate allele, so that it is informative!
			if (vd.getGenotypeCountReference()>0 || vd.getGenotypeCountAlternate()>0)
				return (vd);
		}		
		return null;
	}

	/**
	 * Computes a missing data penalty for each pileup/variant.
	 * This penality is based on the genotype clases observed for this SNP across all samples, and is a weighted mix of
	 * L(D|HOM_REF)*fraction HOM_REF + L(D|HET)*fraction HET + L(D|HOM_VAR)*fraction HOM_VAR
	 * The value returned is NOT in log space.
	 * @param allProbs
	 * @param genotypeMatrix
	 * @param fixedGenotypeErrorRate
	 * @return
	 */
	public Map<SampleGenotypeProbabilities, Double> getMissingDataPenalities(final List<SampleGenotypeProbabilities> allProbs, final GenotypeMatrix genotypeMatrix, final Double fixedGenotypeErrorRate, final Double maximumObservationProbability) {
		Map<SampleGenotypeProbabilities, Double> result = new HashMap<>();
		for (SampleGenotypeProbabilities p: allProbs) {
			Interval i = p.getSNPInterval();
			Double minorAlleleFrequency = CellAssignmentUtils.getNullableValue(this.variantMinorAlleleFrequency, i);						
			double [] genotypeFreqs = genotypeMatrix.getGenotypeFrequencies(i);
			char refAllele = StringUtil.byteToChar(genotypeMatrix.getRefBase(i));
			char altAllele = StringUtil.byteToChar(genotypeMatrix.getAltBase(i));
			GenotypeType [] genos = {GenotypeType.HOM_REF, GenotypeType.HET, GenotypeType.HOM_VAR};

			List<Double> freqs=new ArrayList<>();
			for (double genotypeFreq : genotypeFreqs)
				freqs.add(genotypeFreq);
			
			double penalty = p.getLogLikelihoodMissingData(refAllele, altAllele, Arrays.asList(genos), freqs, fixedGenotypeErrorRate, null, maximumObservationProbability, minorAlleleFrequency, this.contamination);
			// go back to normal space instead of log space to be consistent with other data.
			penalty=Math.pow(10,penalty);
			result.put(p, penalty);
		}
		return (result);
	}



	public String getCell() {
		return cell;
	}
}
