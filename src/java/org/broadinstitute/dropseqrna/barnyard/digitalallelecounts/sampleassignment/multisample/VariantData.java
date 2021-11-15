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

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.GenotypeType;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.LikelihoodUtils;

public class VariantData {

	private final Interval snpInterval;
	private final byte refAllele;
	private final byte altAllele;
	private final List<Byte> bases;
	private final List<Byte> qualities;
	private final List<GenotypeType> genotypes;
	private final Double missingDataPenality;
	private final Double maximumObservationProbability;
	private final double [] [] allLikes;
	private final Double minorAlleleFrequency;
	private final Double contamination;
	
	public VariantData(final Interval snpInterval, final byte refAllele, final byte altAllele,
			final GenotypeType s1, final GenotypeType s2,
			final List<Byte> bases, final List<Byte> qualities, final Double missingDataPenality, final Double maximumObservationProbability,
			final Double contamination, final Double minorAlleleFrequency) {
		
		this.snpInterval=snpInterval;
		this.refAllele = refAllele;
		this.altAllele = altAllele;
		this.bases = bases;
		this.qualities = qualities;
		this.genotypes = Arrays.asList(s1, s2);
		this.missingDataPenality=missingDataPenality;
		this.maximumObservationProbability=maximumObservationProbability;
		this.minorAlleleFrequency=minorAlleleFrequency;
		this.contamination=contamination;
		this.allLikes= precomputeLikelihoods();
	}

	public VariantData(final Interval snpInterval, final byte refAllele, final byte altAllele,
			final GenotypeType s1, final GenotypeType s2,
			final List<Byte> bases, final List<Byte> qualities) {

		this(snpInterval, refAllele, altAllele, s1, s2, bases, qualities, null, null, null, null);
	}

	
	/**
	 * A convenience constructor that uses more human-readable objects.
	 * 
	 * @param snpInterval The genomic coordinates of the SNP
	 * @param refAllele The reference allele of the observed SNP in the pileup
	 * @param altAllele The alternate allele of the observed SNP in the pileup
	 * @param s1 Genotype of sample 1
	 * @param s2 Genotype of sample 2
	 * @param bases The bases observed in the sequencing pileup
	 * @param qualities The qualities of the bases observed in the sequencing pileup.
	 */
	public VariantData (final Interval snpInterval, final char refAllele, final char altAllele,
			final GenotypeType s1, final GenotypeType s2,
			final char [] bases, final int [] qualities, final Double missingDataPenality, final Double maximumObservationProbability, 
			final Double contamination, final Double minorAlleleFrequency) {

		if (bases.length!=qualities.length)
			throw new IllegalArgumentException("Number of bases and qualities are different!");
		byte ref = StringUtil.charToByte(refAllele);
		byte alt = StringUtil.charToByte(altAllele);
				
		List<Byte> b = new ArrayList<>();
		for (char base: bases)
			b.add(StringUtil.charToByte(base));

		List<Byte> q = new ArrayList<>();
		for (int qual: qualities)
			q.add((byte)qual);

		this.snpInterval=snpInterval;
		this.refAllele=ref;
		this.altAllele=alt;
		this.bases=b;
		this.qualities=q;
		this.genotypes = Arrays.asList(s1, s2);
		this.missingDataPenality=missingDataPenality;
		this.maximumObservationProbability=maximumObservationProbability;
		this.minorAlleleFrequency=minorAlleleFrequency;
		this.contamination=contamination;		
		this.allLikes= precomputeLikelihoods();
	}
	
	/**
	 * Since the genotype states, ref/alt allele, and bases/qualities are fixed, compute the likelihood of each genotype for each UMI pileup base observed.
	 * 
	 * @return A 2d matrix.  The first dimension is 1 base of the UMI pileup.  The second dimension contains the likelihoods for the genotypes given the UMI base/quality.
	 */
	private double [] [] precomputeLikelihoods () {
		double [] [] allLikes = new double [bases.size()] [genotypes.size()];
		for (int i=0; i<bases.size(); i++) {
			//TODO: add in contamination calculations
			allLikes[i] = LikelihoodUtils.getInstance().getLikelihoodManyObservations(this.refAllele, this.altAllele, this.genotypes, this.bases.get(i), this.qualities.get(i),
					this.missingDataPenality, null, this.maximumObservationProbability, this.minorAlleleFrequency, this.contamination);
		}
		return (allLikes);
		
	}

	
	public VariantData (final Interval snpInterval, final char refAllele, final char altAllele,
			final GenotypeType s1, final GenotypeType s2,
			final char [] bases, final int [] qualities) {
		this(snpInterval, refAllele, altAllele, s1, s2, bases, qualities, null, null, null, null);
	}
	
	public Interval getSNPInterval() {
		return this.snpInterval;
	}


	/**
	 * Is genotype one different from genotype two, and is neither of the genotypes a no call?
	 * @return
	 */
	public boolean isInformative() {
		if (this.getGenotypeOne()==GenotypeType.NO_CALL || this.getGenotypeTwo()==GenotypeType.NO_CALL) return false;
		if (this.getGenotypeOne()==this.getGenotypeTwo()) return false;
		return true;
	}
	
	/**
	 * Is genotype one different from genotype two, and is neither of the genotypes a no call?
	 * Additionally, the genotype for the sample should be homozygous. 
	 * @param genotypeIndex 
	 * @return True if this 
	 */
	public boolean isInformativeHomozygous(final int genotypeIndex) {
		if (!isInformative()) return false;
		// Site is informative, is it homozygous?
		GenotypeType gt = this.genotypes.get(genotypeIndex);
		switch (gt) {
		case HOM_REF: return true;
		case HOM_VAR: return true;
		default: return false;
		}
	}
	
	public int getNumImpossibleAllelesGenotypeOne() {
		return getNumImpossibleAlleles(0);
	}

	public int getNumImpossibleAllelesGenotypeTwo() {
		return getNumImpossibleAlleles(1);
	}

	int getNumImpossibleAlleles(final int genotypeIndex) {
		GenotypeType gt = this.genotypes.get(genotypeIndex);
		int count=0;
		switch (gt) {
			case HOM_REF:
				count = getCountAllele(this.altAllele);
				break;
			case HET:
				count=0;
				break;
			case HOM_VAR:
				count = getCountAllele(this.refAllele);
				break;
		}
		return count;
	}

	private int getCountAllele (final byte allele) {
		int count=0;
		for (byte b: this.bases)
			if (b==allele)
				count++;
		return (count);
	}

	public GenotypeType getGenotypeOne () {
		return this.genotypes.get(0);
	}

	public GenotypeType getGenotypeTwo () {
		return this.genotypes.get(1);
	}

	public int getGenotypeCountReference() {
		return getCountAllele(this.refAllele);
	}

	public int getGenotypeCountAlternate() {
		return getCountAllele(this.altAllele);
	}

	/**
	 * How many observations are there of any allele?
	 * @return
	 */
	public int getNumUMIs() {
		int result = getGenotypeCountReference() + getGenotypeCountAlternate();
		return result;
	}

	public double getLogLikelihood(final double mixture) {
		List<Double> mixtures = Arrays.asList(mixture, 1-mixture);
//		double perSNP = LikelihoodUtils.getInstance()
//				.getLogLikelihoodMixedModel(refAllele, altAllele, genotypes,
//						mixtures, bases, qualities, missingDataPenality, null, maximumObservationProbability);
		double perSNP = LikelihoodUtils.getInstance().getLogLikelihoodMixedModel(this.allLikes, mixtures);
//		if (Math.abs(perSNP-perSNP2)>0.001)
//			System.out.println("STOP");
		if (Double.isNaN(perSNP))
			throw new IllegalStateException("Likelihood result can't be NaN");
		return perSNP;
	}

	@Override
	public String toString () {
		StringBuilder b = new StringBuilder();
		b.append(snpInterval.toString()+" ");
		b.append("Ref [" +StringUtil.byteToChar(this.refAllele)+"] ");
		b.append("Alt [" +StringUtil.byteToChar(this.altAllele)+"] ");
		b.append("Genotype1 ["+ this.genotypes.get(0).toString()+"] ");
		b.append("Genotype2 ["+ this.genotypes.get(1).toString()+"] ");
		b.append("Bases [" +StringUtil.bytesToString(convert(this.bases))+"] ");
		b.append("Quals [");
		for (int i=0; i<this.qualities.size(); i++) {
			b.append(this.qualities.get(i));
			if (i<this.qualities.size()-1)
				b.append(",");
		}
		b.append("]");
		return b.toString();
	}

	private byte [] convert (final List<Byte> x) {
		byte[] base =new byte[x.size()];
		for (int i=0; i<x.size(); i++)
			base[i]=x.get(i);
		return base;
	}

}
