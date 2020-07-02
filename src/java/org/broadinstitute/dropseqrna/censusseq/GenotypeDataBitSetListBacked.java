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
package org.broadinstitute.dropseqrna.censusseq;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import htsjdk.variant.variantcontext.GenotypeType;

public class GenotypeDataBitSetListBacked implements GenotypeDataI {

	// one bitset per sample.
	private List<BitSet> data;

	private final int BITS_PER_DATUM=2;

	public GenotypeDataBitSetListBacked (final int numSamples) {
		// initialize a bitset for each sample.
		data=new ArrayList<>(numSamples);
		for (int i=0; i<numSamples; i++)
			data.add(new BitSet());
	}

	@Override
	public void add(final int variantIndex, final int sampleIndex, final int altAlleleCount) {
		int modifier=encode(altAlleleCount);
		if (modifier==-1) return;

		int idx= getBitsStartLocation(variantIndex);
		BitSet sampleBitSet = data.get(sampleIndex);
		if (modifier==-2) {
			sampleBitSet.set(idx);
			sampleBitSet.set(idx+1);
			return;
		}
		sampleBitSet.set(idx+modifier);
	}
	
	public void add(final int variantIndex, final int sampleIndex, final GenotypeType genotype) {
		int altAlleleCount = genotypeToAlleleCount(genotype);
		add(variantIndex, sampleIndex, altAlleleCount);
	}

	@Override
	public int get(final int variantIndex, final int sampleIndex) {
		int idx= getBitsStartLocation(variantIndex);
		BitSet sampleBitSet = data.get(sampleIndex);
		int result = decode(sampleBitSet, idx, idx+(BITS_PER_DATUM-1));
		return result;
	}
	
	public GenotypeType getAsGenotypeType (final int variantIndex, final int sampleIndex) {
		int altAlleleCount = get (variantIndex, sampleIndex);
		return (altAlleleCountToGenotypeType(altAlleleCount));
	}

	/**
	 * Store 2 bits per result.
	 * @param variantIndex
	 * @param sampleIndex
	 * @return
	 */
	private int getBitsStartLocation (final int variantIndex) {
		int idx1 = (variantIndex)*BITS_PER_DATUM;
		return (idx1);
	}

	/**
	 * If both bits are set, data is missing (-1).
	 * If neither bit is set, number of alt alleles = 0.
	 * If first bit is set, number of alt alleles =1;
	 * If second bit is set, number of alt alleles =2;
	 * @param startIndex
	 * @param endIndex
	 * @return
	 */
	private int decode (final BitSet sampleBitSet, final int startIndex, final int endIndex) {
		boolean f1 = sampleBitSet.get(startIndex);
		boolean f2 = sampleBitSet.get(endIndex);
		if (f1 && f2) return -1;
		if (f1) return 1;
		if (f2) return 2;
		return 0;
	}
		
	/**
	 * See {@link}decode
	 * returns -1 if nothing should be set, 0 to set the first bit, or 1 to set the second.
	 * This corresponds to what bit should be toggled from the starting index.
	 * For 0 alt alleles, if your starting index is 0 you receive a -1 and do nothing.
	 * For 1 alt allele, if your starting index is 0 you receive a 0, and set the first bit in the bitset.
	 * For 2 alt alleles, if your starting index is 0 you receive a 1, and set the second bit 1 in the bitset.
	 * For missing data (-1), set both bits.
	 * @param altAlleleCount The number of alternate alleles for this sample/variant.
	 */
	private int encode (final int altAlleleCount) {
		if (altAlleleCount==-1) return -2;
		return altAlleleCount-1;
	}
	
	/**
	 * Convert a GenotypeType to a count of alternate alleles.
	 * HOM_REF = 0 alt alleles
	 * HET = 1 alt alleless
	 * HOM_VAR = 2 alt alleles
	 * Anything else = missing.
	 * @param g the genotype
	 * @return the number of alternate alleles.
	 */
	private int genotypeToAlleleCount (final GenotypeType g) {
		switch (g) {
		case HOM_REF:return (0);
		case HET: return (1);
		case HOM_VAR: return (2);
		default: return (-1);
		}					
	}
	
	/*
	 * Convert a count of alt alleles to a GenotypeType.
	 * All data that wasn't HOM_REF/HET/HOM_VAR is encoded as missing.
	 * @param altAlleleCount the number of alternate alleles in the genotype
	 * @return the Genotype for the number of alternate alleles.
	 */
	private GenotypeType altAlleleCountToGenotypeType (final int altAlleleCount) {
		switch (altAlleleCount) {
		case 0: return (GenotypeType.HOM_REF);
		case 1: return (GenotypeType.HET);
		case 2: return (GenotypeType.HOM_VAR);
		default: return (GenotypeType.NO_CALL);
	}
	}
	
	
}
