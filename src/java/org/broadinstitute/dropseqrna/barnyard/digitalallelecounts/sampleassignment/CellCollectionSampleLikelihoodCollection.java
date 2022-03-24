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

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.dropseqrna.utils.statistics.FDR;

import com.google.common.collect.ImmutableMap;

import htsjdk.samtools.util.Interval;
import htsjdk.variant.variantcontext.GenotypeType;
import picard.util.TabbedInputParser;

/**
 * A collection of CellSampleLikelihoodCollection objects for a single variant position.
 * @author nemesh
 *
 */
public class CellCollectionSampleLikelihoodCollection {
	// the key is the cell ID, the value is the collection of likelihoods to samples.
	private final Map<String, CellSampleLikelihoodCollection> map;
	
	// this may be null.
	private final Double fixedGenotypeErrorRate;
	
	// this may be null.
	final Double maximumObservationProbability;
	
	// Map of the donor name to the position of that donor name in the list.  Shared by all CellSampleListlihoodCollection objects.
	private final ImmutableMap<String, Integer> sampleIndexMap;
	private final  Map<String,Double> cellContaminationMap;
	private final Map<Interval, Double> variantMinorAlleleFrequency;	
	private static final GenotypeType [] genotypeModels = {GenotypeType.HOM_REF, GenotypeType.HET, GenotypeType.HOM_VAR};

	
	public CellCollectionSampleLikelihoodCollection(List<String> sampleList) {
		this(null, (Double) null, sampleList);
	}

	/**
	 * Build a collection that holds one or more cells, and the likelihood of the cell's genotype data having come from one or more donors.
	 * @param fixedGenotypeErrorRate Instead of using the base qualities to compute likelihoods, used a fixed error rate instead.
	 * if set to null, ignore this parameter and use the base error rates.
	 * @param maximumObservationProbability Cap the maximum error penalty at this number.
	 */
	public CellCollectionSampleLikelihoodCollection(final Double fixedGenotypeErrorRate, final Double maximumObservationProbability, List<String> sampleList) {
		map = new HashMap<>();
		this.fixedGenotypeErrorRate=fixedGenotypeErrorRate;
		this.maximumObservationProbability=maximumObservationProbability;
		this.variantMinorAlleleFrequency=null;
		this.cellContaminationMap=null;
		// build the map of the sample list to the positions, retaining order.
		LinkedHashMap<String, Integer> sampleIndexMap = new LinkedHashMap<>();
		for (int i=0; i<sampleList.size(); i++) {
			sampleIndexMap.put(sampleList.get(i), i);
		}
		this.sampleIndexMap=ImmutableMap.copyOf(sampleIndexMap);		
	}
	
	/**
	 * Build a collection that holds one or more cells, and the likelihood of the cell's genotype data having come from one or more donors.
	 */
	public CellCollectionSampleLikelihoodCollection(final Map<String,Double> cellContaminationMap, Map<Interval, Double> variantMinorAlleleFrequency, List<String> sampleList) {
		map = new HashMap<>();
		this.fixedGenotypeErrorRate=null;
		this.maximumObservationProbability=null;
		this.variantMinorAlleleFrequency=variantMinorAlleleFrequency;
		this.cellContaminationMap=cellContaminationMap;
		// build the map of the sample list to the positions, retaining order.
		LinkedHashMap<String, Integer> sampleIndexMap = new LinkedHashMap<>();
		for (int i=0; i<sampleList.size(); i++) {
			sampleIndexMap.put(sampleList.get(i), i);
		}
		this.sampleIndexMap=ImmutableMap.copyOf(sampleIndexMap);				
	}
	
	/**
	 * If a contamination score was set for this cell, return it, else return null
	 * @param cellBarcode The cell barcode to query.
	 * @return The contamination score, or null
	 */
	public Double getContamination (String cellBarcode) {
		if (this.cellContaminationMap==null) return (null);
		return (this.cellContaminationMap.get(cellBarcode));
	}
	
	public boolean hasContamination () {
		if (this.cellContaminationMap==null) return (false);
		return !this.cellContaminationMap.isEmpty();
	}
	
	
	// Constructors purely for testing Object.  
	CellSampleLikelihoodCollection buildCellSampleLikelihoodCollection (final String cellBarcode) {
		return (new CellSampleLikelihoodCollection(cellBarcode, 0, 0, 0d, this.sampleIndexMap));
	}
	
	CellSampleLikelihoodCollection buildCellSampleLikelihoodCollection (final String cellBarcode, final int numSNPs, final int numUMIs, final double globalLikelihoodScore) {
		return (new CellSampleLikelihoodCollection(cellBarcode, numSNPs, numUMIs, globalLikelihoodScore, this.sampleIndexMap));
	}
	
	/**
	 * Merges a CellCollectionSampleLikelihoodCollection into this obect. 
	 * @param other Another CellCollectionSampleLikelihoodCollection object.
	 */
	/*  This functionality is disabled as it is unused.
	public void add (final CellCollectionSampleLikelihoodCollection other) {
		if (this.fixedGenotypeErrorRate!=other.fixedGenotypeErrorRate)
			throw new IllegalArgumentException("Added a collection in with a different fixed error rate.");
		for (String cellBarcode: other.getCellBarcodes()) {
			CellSampleLikelihoodCollection c = this.getLikelihoodCollection(cellBarcode);
			if (c==null)
				c = other.getLikelihoodCollection(cellBarcode);
			else
				c.add(other.getLikelihoodCollection(cellBarcode));
			this.map.put(cellBarcode, c);
		}
	}
	*/
	
	/**
	 * Add a CellSampleLikelihoodCollection to this collection.
	 */
	void add (final CellSampleLikelihoodCollection data) {
		this.map.put(data.getCellBarcode(), data);
	}		

	/**
	 * Get the set of samples from the VCF represented anywhere in the data.
	 * @return
	 */
		
	/**
	 * Update the number of SNPs and UMIs observed for the cell captured by the pileup.
	 * Call this for each cell/snp pairing.
	 */
	public void updateNumObservations(final SampleGenotypeProbabilities p) {
		String cellBarcode = p.getCell();
		CellSampleLikelihoodCollection c = this.map.get(cellBarcode);
		if (c==null)
			c=new CellSampleLikelihoodCollection(cellBarcode, 0, 0, 0, this.sampleIndexMap);
		int numUMIs = p.getBackingPileups().size();
		c.incrementNumSNPs(1);
		c.incrementNumUMIs(numUMIs);
		// this should be part of the null check.
		this.map.put(cellBarcode, c);
	}
						
		
	
	/**
	 * For a pileup of reads on a cell/snp and a sample with two allele calls (ie: a genotype)
	 * update the likelihood that the cell should be assigned to the sample.
	 * 
	 * 
	 * @param p A pileup of observations at a SNP for a single cell.
	 * @param sample The sample being updated
	 * @param alleleOne The first allele of the SNP for the sample
	 * @param alleleTwo The second allele of the SNP for the sample
	 * @param genotypeProbability The probability of the genotype.  Can be set to null to not consider.
	 * @param referenceAllele The reference allele for this variant in the population. (Optional, set null to ignore).
	 */
	public void updatelikelihoods (final SampleGenotypeProbabilities p, final String sample, final char alleleOne, final char alleleTwo, final Double genotypeProbability, final Character referenceAllele) {
		Double maf = CellAssignmentUtils.getNullableValue(this.variantMinorAlleleFrequency, p.getSNPInterval());
		Double contamination = CellAssignmentUtils.getNullableValue(this.cellContaminationMap, p.getCell());
		double logLikelihood = p.getLogLikelihood(alleleOne, alleleTwo, this.fixedGenotypeErrorRate, genotypeProbability, this.maximumObservationProbability, referenceAllele, maf, contamination);
		updatelikelihoods(p, sample, logLikelihood);
	}
		
	private void updatelikelihoods(final SampleGenotypeProbabilities p, final String sample, double logLikelihood) {
		String cellBarcode = p.getCell();		
		CellSampleLikelihoodCollection c = this.map.get(cellBarcode);
		if (c==null) {
			c=new CellSampleLikelihoodCollection(cellBarcode, 0, 0, 0, this.sampleIndexMap);
			this.map.put(cellBarcode, c);
		}
		c.addLikelihood(sample, logLikelihood);
	}
	
	
	/**
	 * Update the likelihoods for this SNP for samples that do not have this genotype set.
	 * To do this, create a blended model of the population of known genotypes to calculate a new genotype probability.
	 * Apply this likelihood to each unknown sample.
	 * @param p the pileup of read observations to process.
	 * @param refAllele the reference allele for this variant
	 * @param altAllele the alternate allele for this variant
	 * @param homRefCount how many homozygous reference donors were observed for this variant
	 * @param hetCount how many heterozygous donors were observed for this variant
	 * @param altCount how many homozygous alternate donors were observed for this variant
	 * @param samples the list of missing samples.
	 * @param genotypeProbability The probability of the genotype, between 0 and 1 inclusive.  Can be set to null to ignore.
	 */
//	
//	public void updateMissingLikelihoods (final SampleGenotypeProbabilities p, final char refAllele, final char altAllele, final int homRefCount, final int hetCount, final int homVarCount, 
//			final Collection<String> samples, final Double genotypeProbability) {
//		double logLikelihood = getMissingLikelihood(p, refAllele, altAllele, homRefCount, hetCount, homVarCount, genotypeProbability);
//		applyMissingLikelihoodToDonors(p, logLikelihood, samples);				
//	}
//		
	
	/**
	 * Update the likelihoods for this SNP for samples that do not have this genotype set.
	 * To do this, create a blended model of the population of known genotypes to calculate a new genotype probability.
	 * Apply this likelihood to each unknown sample.
	 * @param p the pileup of read observations to process.
	 * @param logLikelihood the missing data penalty score to apply to all samples.
	 * @param samples the list of missing samples.
	 * 
	 */
	public void applyMissingLikelihoodToDonors (final SampleGenotypeProbabilities p, double logLikelihood, final Collection<String> samples) {
		String cellBarcode = p.getCell();
		CellSampleLikelihoodCollection c = this.map.get(cellBarcode);
		if (c==null) {
			c=new CellSampleLikelihoodCollection(cellBarcode, 0, 0, 0, this.sampleIndexMap);
			this.map.put(cellBarcode, c);
		}
		for (String sample: samples)
			c.addLikelihood(sample, logLikelihood);
	}
	
	/**
	 * For each cell/SNP, increase the global likelihood penalty score.
	 * @param p the pileup of read observations to process.
	 * @param logLikelihood the missing data penalty score
	 */
	public void incrementGlobalPenaltyAllDonors(final SampleGenotypeProbabilities p, double logLikelihood) {
		String cellBarcode = p.getCell();
		CellSampleLikelihoodCollection c = this.map.get(cellBarcode);
		if (c==null) {
			c=new CellSampleLikelihoodCollection(cellBarcode, 0, 0, 0, this.sampleIndexMap);
			this.map.put(cellBarcode, c);
		}
		c.addGlobalPenalty(logLikelihood);		
	}
	
	/**
	 * Create a blended model of the population of known genotypes to calculate a new genotype probability.
	 * @param p the pileup of read observations to process.
	 * @param refAllele The first allele of the SNP for the sample
	 * @param altAllele The second allele of the SNP for the sample
	 * @param homRefCount how many homozygous reference donors were observed for this variant
	 * @param hetCount how many heterozygous donors were observed for this variant
	 * @param genotypeProbability The probability of the genotype, between 0 and 1 inclusive.  Can be set to null to ignore.
	 * @return The log likelihood of the data for the blended model.
	 */
	public double getMissingLikelihood(final SampleGenotypeProbabilities p, final char refAllele, final char altAllele, final int homRefCount, final int hetCount, final int homVarCount, 
			final Double genotypeProbability) {
		
		Double minorAlleleFrequency = CellAssignmentUtils.getNullableValue(this.variantMinorAlleleFrequency, p.getSNPInterval());
		Double contamination = CellAssignmentUtils.getNullableValue(this.cellContaminationMap, p.getCell());
		List<GenotypeType> genotypes  = Arrays.asList(genotypeModels);
		List<Double> mixture = Arrays.asList((double) homRefCount, (double) hetCount, (double) homVarCount);
		double logLikelihood = p.getLogLikelihoodMissingData(refAllele, altAllele, genotypes, mixture, this.fixedGenotypeErrorRate, genotypeProbability, this.maximumObservationProbability, minorAlleleFrequency, contamination);
		return logLikelihood;
	}

	/**
	 * Get a set of cell barcodes contained in this dataset.
	 */
	public Set<String> getCellBarcodes () {
		return this.map.keySet();
	}

    /**
     * Get sample list contained in this dataset.
     * @return List of samples
     */
    public List<String> getSamples() {
        return new ArrayList<>(this.sampleIndexMap.keySet());
    }

	public CellSampleLikelihoodCollection getLikelihoodCollection (final String cellBarcode) {
		return this.map.get(cellBarcode);
	}
		
	/**
	 * After all data has been added (so you have the best donor per cell), correct the collection of donors pvalues by FDR.
	 * Return a map of each cell ID to the FDR corrected pvalue.
	 */
	public Map<String, Double> getFDRCorrectedPvalues () {
		Iterator<String>cellBCs=this.map.keySet().iterator();
		double [] pvalues = new double [this.map.size()];

		for (int i=0; i<pvalues.length; i++) {
			String cellBC=cellBCs.next();
			CellSampleLikelihoodCollection c = this.map.get(cellBC);
			pvalues[i] = c.getBestSampleAssignment().getOneMinusPvalue();
		}
		double [] fdrCorrectedPvalues = FDR.bhAdjustment(pvalues);

		cellBCs=this.map.keySet().iterator();
		Map<String, Double> result = new HashMap<>();
		for (int i=0; i<pvalues.length; i++) {
			String cellBC=cellBCs.next();
			result.put(cellBC, fdrCorrectedPvalues[i]);
		}
		return result;

	}

	public static CellCollectionSampleLikelihoodCollection parseFile(final File data) {
		
		TabbedInputParser parser = new TabbedInputParser(false, data);
		Iterator<String [] > iter = parser.iterator();
		String [] header =iter.next();
		// the expectation is that the first sample starts after either the bestSample or the trueSample column.
		int indexSampleStart=getIndexOfFirstSample(header);

		// test if there is a global likelihood score in the header.
		int palIdx=Arrays.asList(header).indexOf("population_average_likelihood");
		
		String [] samples = Arrays.copyOfRange(header, indexSampleStart,  header.length);		
		CellCollectionSampleLikelihoodCollection result = new CellCollectionSampleLikelihoodCollection(Arrays.asList(samples));
		// parse each line into a CellSampleLikelihoodCollection object.
		while (iter.hasNext()) {
			String [] line = iter.next();
			double pal = 0;
			if (palIdx!=-1) pal=Double.parseDouble(line[palIdx]);
			CellSampleLikelihoodCollection o = result.buildCellSampleLikelihoodCollection(line[0],Integer.parseInt(line[1]), Integer.parseInt(line[2]), pal);
			String [] likes = Arrays.copyOfRange(line, indexSampleStart,  line.length);
			for (int i=0; i<likes.length; i++) {
				String sample = samples[i];
				double logLikelihood = Double.parseDouble(likes[i]);
				o.addLikelihood(sample, logLikelihood);
			}
			result.add(o);

		}

		parser.close();

		return result;
	}

	private static int getIndexOfFirstSample (final String [] header) {
		// any one of these can be the last column before the samples, depending on file versioning and if it's positive control data.
		List<String> lastColumn = Arrays.asList("bestSample", "trueSample", "median_likelihood", "population_average_likelihood");
		
		int index =-1;
		for (int i=0; i<header.length; i++) {
			String column=header[i];
			if (lastColumn.contains(column)) {
				if (i>index) index=i;
			}
		}			
		return index+1;
	}
	
	
}
