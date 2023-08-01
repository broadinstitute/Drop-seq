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
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPInfoCollection;
import org.broadinstitute.dropseqrna.censusseq.GenotypeDataBitSetListBacked;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.vcftools.filters.GenotypeGQFilter;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContext.Type;

/**
 * Holds genotypes at intervals across samples.
 * @author nemesh
 *
 */
public class GenotypeMatrix {

	private static final Log log = Log.getInstance(GenotypeMatrix.class);

	private final GenotypeDataBitSetListBacked matrix;
	private Map<Interval, Integer> rows;
	private Map<String, Integer> columns;
	// calculate this on the fly and cache.
	private Map<Interval, double []> genotypeAlleleFractionMap;
	private final List<Byte> refAllelles;
	private final List<Byte> altAllelles;
	// this is similar to the SNPInfoCollection.
	private Map<Interval, Double> averageGQ;
	
	public GenotypeMatrix(final GenotypeDataBitSetListBacked matrix, final Map<Interval, Integer> rows, final Map<String, Integer> columns, final List<Byte> refAllelles, final List<Byte> altAllelles) {
		this.matrix=matrix;
		this.rows=rows;
		this.columns=columns;
		this.refAllelles=refAllelles;
		this.altAllelles=altAllelles;
		this.genotypeAlleleFractionMap=new HashMap<>();
		this.averageGQ = new HashMap<>();
		// cache genotype frequencies.
		rows.keySet().forEach(x -> getGenotypeFrequencies(x));
	}
		
	/**
	 * Build a genotype matrix based on the variant context iterator.
	 * Each variant must have at least one sample called GenotypeType.HOM_REF, GenotypeType.HET, or GenotypeType.HOM_VAR
	 * @param iter a variant context iterator, pre-filter this as much as you'd like
	 * @param genotypeQuality The minimum genotype quality to set.  Set this to null to not filter on genotype quality.
	 */
	public GenotypeMatrix(final Iterator<VariantContext> iter, final Integer genotypeQuality, final Collection<String> donors) {
		rows = new HashMap<>();
		columns = new HashMap<>();
		Set<String> sampleSet = new HashSet<>(donors);

		// Map<Integer, Map<Integer, GenotypeType>> tempMatrix = new HashMap<>();
		this.matrix = new GenotypeDataBitSetListBacked(donors.size());
		this.genotypeAlleleFractionMap=new HashMap<>();
		this.averageGQ = new HashMap<>();
		
		this.refAllelles = new ArrayList<>();
		this.altAllelles = new ArrayList<>();
		int count=0;
		int maxRows=-1;

		log.info("Adding VCF genotypes to Genotype Matrix in memory");
		while (iter.hasNext()) {
			
			count++;
			if (count%100000==0) log.info("Added " + count + " variant records");
			VariantContext vc = iter.next();
			Interval i = new Interval(vc.getContig(), vc.getStart(), vc.getEnd());
			
			// calculate and store the average site genotype quality for later variant filtering.
			double gqAverage = vc.getGenotypes().stream().mapToInt(x -> x.getGQ()).average().orElse(-1);
			this.averageGQ.put(i, gqAverage);
			
			GenotypeGQFilter gf = new GenotypeGQFilter(vc.getGenotypes().iterator(), genotypeQuality);			
			List<String> samples = vc.getSampleNamesOrderedByName();
			byte refBase = vc.getReference().getBases()[0];
			byte altBase = getAltBase(vc);
			refAllelles.add(refBase);
			altAllelles.add(altBase);
			for (String sample: samples) {
				if (!sampleSet.contains(sample))
					continue;
				Genotype g = vc.getGenotype(sample);
				GenotypeType t = g.getType();
				if (t==null) // this should never happen.
					log.info("STOP");
				// put a no call in this position if the genotype should be filtered out.
				if (gf.filterOut(g))
					t=GenotypeType.NO_CALL;
				// get row index.
				Integer rowIndex= rows.get(i);
				if (rowIndex==null) {
					// get the max row and increment.
					maxRows++;
					// Holy crap this gets expensive.  Just track with an incrementor instead of counting objects, which scales awfully.
					//rowIndex=getMaxIndex(rows);
					rowIndex=maxRows;
					rows.put(i, rowIndex);
					if (rowIndex!=maxRows)
						log.info("STOP");

				}
				// get column index.
				Integer colIndex = columns.get(sample);
				if (colIndex==null) {
					colIndex=getMaxIndex(columns);
					columns.put(sample, colIndex);
				}
				//add the data
				matrix.add(rowIndex, colIndex, t);
				
			}
			// cache all the genotye frequencies.
			getGenotypeFrequencies(i);
			gf.close();
		}
				
		log.info("Added " + count + " variant records");
		log.info("Finished adding VCF genotypes to Genotype Matrix in memory");
						
	}

	
	
	
	public Set<Interval> getSNPIntervals () {
		return rows.keySet();
	}
	
	public Byte getRefBase (final Interval i) {
		Integer rowIdx= this.rows.get(i);
		if (rowIdx==null) return null;
		return this.refAllelles.get(rowIdx);
	}

	public Byte getAltBase (final Interval i) {
		Integer rowIdx= this.rows.get(i);
		if (rowIdx==null) return null;
		return this.altAllelles.get(rowIdx);
	}

	private Byte getAltBase (final VariantContext vc) {
		// short circuit.
		if (vc.getType()==Type.NO_VARIATION)
		return StringUtil.charToByte('N');
		// otherwise there's variation, continue.
		Allele alt = vc.getAltAlleleWithHighestAlleleCount();
		byte altBase=StringUtil.charToByte('N');

		if (alt!=null) {
			byte [] altBases = alt.getBases();
			if (altBases.length>0)
				altBase=altBases[0];
		}
		return altBase;
	}

	private <T> int getMaxIndex (final Map<T, Integer> map) {
		if (map.isEmpty()) return 0;
		return Collections.max(map.values())+1;
	}

	public GenotypeType getGenotype (final Interval i, final String sample) {
		Integer rowIndex = this.rows.get(i);
		Integer colIndex = this.columns.get(sample);
		if (rowIndex==null | colIndex==null) return null;
		return this.matrix.getAsGenotypeType(rowIndex, colIndex); 
	}

	/**
	 * Return the fraction of genotypes that are HOM_REF, HET, and HOM_VAR
	 * @param interval
	 * @return An array of doubles representing the fraction of samples HOM_REF, HET, and HOM_VAR in the population of samples in the data set.
	 * returns null if this interval doesn't exist in the data.
	 */
	public double [] getGenotypeFrequencies (final Interval interval) {
		double [] result = this.genotypeAlleleFractionMap.get(interval);
		if (result!=null)
			return result;
		
		if (interval==null) // no interval saved!
			return null;
		
		int [] genotypes = getCountsAltAllele(interval);
		ObjectCounter< Integer> counter = new ObjectCounter<>();
		Arrays.stream(genotypes).forEach(x-> counter.increment(x));
				
		double fracHomRef=(double) counter.getCountForKey(new Integer (0)) / (double) counter.getTotalCount();
		double fracHet=(double) counter.getCountForKey(new Integer (1)) / (double) counter.getTotalCount();
		double fracHomVar=(double) counter.getCountForKey(new Integer (2)) / (double) counter.getTotalCount();
		double [] result2 = {fracHomRef, fracHet, fracHomVar};
		this.genotypeAlleleFractionMap.put(interval, result2);
		return result2;
	}

	
	/**
	 * Gets the genotype states as the count of alternate alleles for a variant.
	 * The return array is in the same order the samples are stored.
	 * @param variantIndex
	 * @return
	 */
	public int [] getCountsAltAllele (final Interval interval) {
		int variantIndex = this.rows.get(interval);
		
		int [] result = new int [this.columns.size()];
		for (int i=0; i<result.length; i++)
			result[i]=this.matrix.get(variantIndex, i);
		return (result);
	}

	/**
	 * Does this data set contain the donor?
	 * @param donor
	 * @return
	 */
	public boolean containsDonor (final String donor) {
		return this.columns.keySet().contains(donor);
	}
	
	/**
	 * The list of intervals stored in the genotype matrix.
	 * @return
	 */
	public Set<Interval> getIntervals () {
		return this.rows.keySet();
	}
	
	public Map<Interval, Double> getAverageGenotypeQuality () {
		return this.averageGQ;
	}

}
