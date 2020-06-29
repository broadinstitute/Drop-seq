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

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.broadinstitute.dropseqrna.utils.ObjectCounter;

import com.google.common.collect.BiMap;
import com.google.common.collect.HashBiMap;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.GenotypeType;
import picard.util.BasicInputParser;

/**
 * Holds the number of ref/alt counts observed across pileups, and the genotype states of all samples.
 * @author nemesh
 *
 */
public class CommonSNPsData {

	private static final Log log = Log.getInstance(CommonSNPsData.class);

	private List <int []> refAltCounts;
	// One BitSet of N*M*2.  Each 2 bits stores an entry.
	// write a bit of storage code that indexes into each row/column and encodes/decodes the two bits -> genotype state.
	// and either inserts them or retrieves them.

	// private List<byte []> genotypeStates;

	private final BiMap<String, Integer> sampleIndexMap;
	private final GenotypeDataI genotypeStates;

	// should missing donors be removed from frequency calculation when their data is missing for a SNP?
	// if true, remove the donor
	// if false, treat the donor as heterozygous.

	private final boolean OMIT_MISSING_DATA=false;

	public CommonSNPsData (final List<String> sampleNames) {
		genotypeStates = new GenotypeDataBitSetListBacked(sampleNames.size());

		refAltCounts = new ArrayList<>();
		sampleIndexMap = HashBiMap.create();

		for (int i=0; i<sampleNames.size(); i++)
			sampleIndexMap.put(sampleNames.get(i), i);
	}

	/**
	 * Add one SNP's worth of data to the data set.
	 * @param sampleNames
	 * @param genotypes
	 * @param refAltCounts
	 */
	public void addSNP (final String [] sampleNames, final int [] altAlleleCounts, final int readsRef, final int readsAlt) {
		int variantIndex = this.refAltCounts.size();
		int [] r = {readsRef, readsAlt};
		this.refAltCounts.add(r);
		for (int i=0; i<sampleNames.length; i++) {
			int sampleIndex = sampleIndexMap.get(sampleNames[i]);
			genotypeStates.add(variantIndex, sampleIndex, altAlleleCounts[i]);
		}
	}

	public ObjectCounter<GenotypeType> getCountGenotypes(final String sample) {
		ObjectCounter<GenotypeType> result = new ObjectCounter<>();
		Integer idxSample = sampleIndexMap.get(sample);
		if (idxSample==null) return null;

		for (int idxVariant=0; idxVariant<getNumVariants(); idxVariant++) {
			int countAltAllele = this.genotypeStates.get(idxVariant, idxSample);
			if (countAltAllele==0) result.increment(GenotypeType.HOM_REF);
			if (countAltAllele==1) result.increment(GenotypeType.HET);
			if (countAltAllele==2) result.increment(GenotypeType.HOM_VAR);
			if (countAltAllele==-1) result.increment(GenotypeType.NO_CALL);
		}
		return result;
	}

	public int getNumAltAlleles (final String sampleName, final int variantIndex) {
		return genotypeStates.get(variantIndex, sampleIndexMap.get(sampleName));
	}

	/**
	 * How many total observations of ref or alt alleles are there across all SNPs?
	 * @return
	 */
	public int getTotalSNPAlleleCounts() {
		int count=0;
		for (int [] d: this.refAltCounts)
			for (int i: d)
				count+=i;
		return (count);
	}


	public List<String> getSampleNames () {
		List<String> result = new ArrayList<>();
		Map<Integer, String> reverseMap = sampleIndexMap.inverse();
		for (int i=0; i<reverseMap.size(); i++)
			result.add(reverseMap.get(i));
		return result;
	}



	public int getNumVariants () {
		return this.refAltCounts.size();
	}

	public int [] getRefAltCounts (final int variantIndex) {
		return this.refAltCounts.get(variantIndex);
	}

	/**
	 * Gets the genotype states as the count of alternate alleles for a variant.
	 * The return array is in the same order the samples are stored.
	 * @param variantIndex
	 * @return
	 */
	public int [] getCountsAltAllele (final int variantIndex) {
		int [] result = new int [sampleIndexMap.size()];
		for (int i=0; i<result.length; i++)
			result[i]=genotypeStates.get(variantIndex, i);
		return (result);
	}

	public double [] getWeighedAlleleFrequencies (final double [] ratios) {
		int numVariants = getNumVariants();
		double [] result = new double [numVariants];
		for (int i=0; i<numVariants; i++)
			result[i]=getWeighedAlleleFrequenciesOneSNP(ratios, i);

		return result;
	}


	public double getWeighedAlleleFrequenciesOneSNP (final double [] ratios, final int variantIndex) {
		if (this.OMIT_MISSING_DATA) return (getWeighedAlleleFrequenciesOneSNPMissingRemoved(ratios, variantIndex));
		return (getWeighedAlleleFrequenciesOneSNPMissingAsHet(ratios, variantIndex));
	}

	/**
	 * Missing data is treated as a heterozygous SNP for this operation.
	 * @param ratios
	 * @param variantIndex
	 * @return
	 */
	public double getWeighedAlleleFrequenciesOneSNPMissingAsHet (final double [] ratios, final int variantIndex) {
		int [] genotypes = getCountsAltAllele(variantIndex);
		double result = 0;
		for (int i=0; i<genotypes.length; i++) {
			// genotypes encoded as -1 are missing
			double gs= genotypes[i];
			if (gs==-1)
				gs =1;
			gs = gs/2;
			result+=ratios[i]*gs;
		}
		if (Double.isNaN(result))
			log.warn("SNP produced NaN allele frequency: ratios [" + ratios.toString() +" ] genotypes [" + this);
		return (result);
	}

	/**
	 * Omit the missing genotypes from the frequency calculation, instead of setting missing data to heterozygous.
	 * @param ratios
	 * @param variantIndex
	 * @return
	 */

	public double getWeighedAlleleFrequenciesOneSNPMissingRemoved (final double [] ratios, final int variantIndex) {
		// if data is missing, need to re-weight the ratios for just this SNP.
		int [] genotypes = getCountsAltAllele(variantIndex);
		return (getWeighedAlleleFrequenciesOneSNPMissingRemoved(ratios, genotypes));

	}


	public double getWeighedAlleleFrequenciesOneSNPMissingRemoved (final double [] ratios, final int [] genotypes) {
		double result = 0;
		int numSkipped=0;
		for (int i=0; i<genotypes.length; i++) {
			double gs= genotypes[i];
			// genotypes encoded as -1 are missing and are skipped.
			if (gs!=-1) {
				gs = gs/2;
				result+=ratios[i]*gs;
			} else
				numSkipped++;
		}
		if (Double.isNaN(result))
			log.warn("SNP produced NaN allele frequency: ratios [" + ratios.toString() +" ] genotypes [" + this);

		if (numSkipped==0)
			return result;

		double fixedResult = getWeighedAlleleFrequenciesOneSNPMissingData(ratios, genotypes);
		return (fixedResult);
	}


	/**
	 * If there's missing genotype states, use this slower method.
	 * Creates the sublist of ratios and genotypes for the non-missing genotype states and calculates allele frequency
	 * based on only non-missing data.
	 * @param ratios
	 * @param variantIndex
	 * @return
	 */
	public double getWeighedAlleleFrequenciesOneSNPMissingData (final double [] ratios, final int [] genotypes) {

		List<Integer> genotypeList = new ArrayList<>(genotypes.length);
		List<Double> ratioList = new ArrayList<>(ratios.length);
		// filter data.
		for (int i=0; i<genotypes.length; i++) {
			int gs= genotypes[i];
			// genotypes encoded as -1 are missing and are skipped.
			if (gs!=-1) {
				genotypeList.add(gs);
				ratioList.add(ratios[i]);
			}
		}
		double[] newRatios = ratioList.stream().mapToDouble(d -> d).toArray();
		// rebalance new ratios to sum to 1 (they will be less than 1 because we dropped at least 1 donor)
		newRatios = OptimizeSampleRatiosGradientFunction.normalizeRatiosToOne(newRatios);
		int [] newGenotypes =  genotypeList.stream().mapToInt(i -> i).toArray();
		return (getWeighedAlleleFrequenciesOneSNPMissingRemoved(newRatios, newGenotypes));
	}

	/**
	 * Parse a pair of files that contains the information needed to optimize donor ratios.
	 *
	 * @param snpReadCountFile A file containing 2 columns (tab sep) with the ref and alt allele counts from the reads.  Has a header.
	 * @param genotypeFile A file containing 1 column per sample, and 1 row per SNP (tab sep).  Has a header with the sample names for the data set.
	 * @return
	 */
	public static CommonSNPsData parseFromFiles (final File snpReadCountFile, final File genotypeFile) {
		IOUtil.assertFileIsReadable(snpReadCountFile);
		IOUtil.assertFileIsReadable(genotypeFile);
		BasicInputParser snpParser = new BasicInputParser(false, snpReadCountFile);
		BasicInputParser genotypeParser = new BasicInputParser(false, genotypeFile);

		// make sure there's data.
		if (!snpParser.hasNext() || !genotypeParser.hasNext()) {
			log.warn("SNP File or Genotype File Empty!");
			snpParser.close();
			genotypeParser.close();
			return null;
		}

		snpParser.next();  // iterate the file parser forward one.
		String [] genotypeHeader= genotypeParser.next(); // get the list of donor names.

		CommonSNPsData result = new CommonSNPsData(Arrays.asList(genotypeHeader));

		while (snpParser.hasNext() & genotypeParser.hasNext()){
			int [] readCounts = parseSNPReadCountLine(snpParser);
			int [] altAlleleCounts = parseSNPReadCountLine(genotypeParser);
			result.addSNP(genotypeHeader, altAlleleCounts, readCounts[0], readCounts[1]);
		}

		snpParser.close();
		genotypeParser.close();
		return (result);
	}



	private static int [] parseSNPReadCountLine (final BasicInputParser parser) {
		if (parser.hasNext()) {
			String [] line =parser.next();
			int [] l = new int [line.length];
			for (int i=0; i<line.length; i++)
				l[i]=Integer.parseInt(line[i]);
			return l;
		}
		return null;
	}


}
