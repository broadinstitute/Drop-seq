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

import htsjdk.samtools.util.Interval;
import org.apache.commons.math3.stat.StatUtils;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.editdistance.MapBarcodesByEditDistance;
import org.broadinstitute.dropseqrna.utils.statistics.BinomialStatistics;

import java.util.*;


// look at this for confidence intervals:
// https://commons.apache.org/proper/commons-math/apidocs/org/apache/commons/math3/stat/interval/IntervalUtils.html

/**
 * For a cell/snp/gene, gather the following results:
 * 1) The number of alleles observed by each UMI, filtered by base quality
 * 2) The number of alleles observed across all UMIs, where each UMI uses the mode allele from step 1
 * 3) Ability to collapse UMIs by edit distance and aggregate those results before calculating step 2
 * 4) ratio of most common vs 2nd most common allele by reads
 * 5) ratio of most common vs 2nd most common allele by UMIs
 * 6) binomial confidence interval (pvalue, lower and upper confidence) for step 4 and step 5.
 *
 * Also store the SNP/gene/cell for the aggregated pileup data.
 * @author nemesh
 *
 */
public class DigitalAlleleCounts {

	private final Interval snpInterval;
	private final String gene;
	private final String cell;
	private final int baseQualityThreshold;
	private Character referenceBase=null;

	private Map<String, ObjectCounter<Character>> umiReadCounts = new HashMap<String, ObjectCounter<Character>>();
	private final MapBarcodesByEditDistance med = new MapBarcodesByEditDistance(false, 1, 0);

	/**
	 * Construct a DAC for a single cell/gene/snp.
	 *
	 * @param snpInterval The SNP interval for this DAC
	 * @param gene The gene for this DAC
	 * @param cell The cell for this DAC
	 * @param baseQualityThreshold Only accept bases that have a quality score of at least baseQualityThreshold.
	 */
	public DigitalAlleleCounts (final Interval snpInterval, final String gene, final String cell, final int baseQualityThreshold) {
		this.snpInterval=snpInterval;
		this.gene=gene;
		this.cell=cell;
		this.baseQualityThreshold=baseQualityThreshold;
	}

	/**
	 * Adds a base pilup which contains the bases and qualities of a SNP on a cell/gene/umi to the DAC collection.
	 * @param p A pileup of bases at a snp/gene/cell/umi.
	 */
	public void addPileup (final SNPUMIBasePileup p) {
		validateIncomingBasePileup(p);
		ObjectCounter<Character> umiBases = convertToObjectCounter(p);
		// a pileup could have all of it's bases filtered out.
		if (umiBases.getSize()>0)
			umiReadCounts.put(p.getMolecularBarcode(), umiBases);
	}

	/**
	 * Add reads to a umi, bypassing base quality threshold.  Adding the same umi string will overwrite data.
	 * This method is far less safe than adding data via a pileup, in that it gives you the freedom to do whatever you'd like
	 * without validating your inputs.
	 * @param umi The umi sequence
	 * @param counts An object counter of base counts for the UMI.
	 */
	public void addReadsForUMI (final String umi, final ObjectCounter<Character> counts) {
		this.umiReadCounts.put(umi, counts);
	}

	/**
	 * Get a list of UMIs that were added to this DAC
	 * @return
	 */
	public Collection <String> umis () {
		return umiReadCounts.keySet();
	}

	/**
	 * Get the allele counts by read for this UMI.
	 * @param umi The name of the UMI
	 * @return
	 */
	public ObjectCounter<Character> getReadsPerUMI (final String umi) {
		return umiReadCounts.get(umi);
	}

	/**
	 * Get the read counts across all UMIs.
	 * This aggregates reads across all possible getReadsPerUMI.
	 * @return
	 */
	public ObjectCounter<Character> getReadCounts() {
		ObjectCounter<Character> result = new ObjectCounter<Character>();
		for (ObjectCounter<Character> o : this.umiReadCounts.values())
			result.increment(o);
		return (result);
	}

	/**
	 * This returns the most common base observed in the read pileup.
	 *
	 * @return
	 */
	public Character getMostCommonBaseByReadCount () {
		return getReadCounts().getKeysOrderedByCount(true).get(0);
	}

	/**
	 * Calculate the confidence interval across all reads on all UMIs in the data set.
	 * This is the way non-UMI aware data sets observe allelic skew.
	 * @param confidenceLevel the desired probability that the true probability of success falls within the returned interval.
	 * @return The binomial confidence interval for the aggregated reads across all UMIs
	 * in this dataset.
	 */
	public BinomialStatistics getReadConfidenceInterval (final double confidenceLevel) {
		ObjectCounter<Character> readCounts = getReadCounts();
		BinomialStatistics result = getBinomialStatistics(readCounts, confidenceLevel);
		return (result);
	}

	public BinomialStatistics getUMIConfidenceInterval (final double confidenceLevel) {
		ObjectCounter<Character> umiCounts = getUMIAlleleCount();
		BinomialStatistics result = getBinomialStatistics(umiCounts, confidenceLevel);
		return (result);
	}


	/**
	 * Calculate the confidence interval for a set of allele counts, using the most common vs second most common alleles.
	 * This uses the Clopper-Pearson approximation, as that's used in R binom.test, which is the standard most commonly used
	 * for this analysis.
	 * @param counts Contains allele counts at either the read or UMI level.
	 * @param confidenceLevel confidenceLevel the desired probability that the true probability of success falls within the returned interval.
	 * @return A ConfidenceInterval object, which contains the upper and lower confidence intervals.
	 */
	public BinomialStatistics getBinomialStatistics (final ObjectCounter<Character> counts, final double confidenceLevel) {
		int [] values = getMostCommonCounts(counts, this.referenceBase);
		int numberOfTrials=values[0]+values[1];
		int numberOfSuccesses=values[0];
		BinomialStatistics bs = new BinomialStatistics(numberOfTrials, numberOfSuccesses, confidenceLevel, 0.5);
		return bs;
	}

	/**
	 * For each UMI, get the mode allele for each UMI.  Sum across all UMIs.
	 * @return
	 */
	public ObjectCounter<Character> getUMIAlleleCount () {
		ObjectCounter<Character> result = new ObjectCounter<Character>();
		for (ObjectCounter<Character> o: this.umiReadCounts.values())
			result.increment(o.getMode());
		return (result);
	}

	/**
	 * This returns the most common base observed in the UMI pileup.
	 *
	 * @return
	 */
	public Character getMostCommonBaseByUMICount () {
		Character result = null;
		try {
			result= getUMIAlleleCount().getKeysOrderedByCount(true).get(0);
		} catch (IndexOutOfBoundsException e) {
			e.printStackTrace();
		}
		return result;
	}

	/**
	 * Get the ratio of the most common base observed by the umi / total number of bases.
	 * @param umi
	 * @return
	 */
	public double getUMIPurity (final String umi) {
		// don't use the reference base, we just want the most common base!
		int [] values = getMostCommonCounts(this.umiReadCounts.get(umi),null);
		int first = values[0];
		int second = values[1];

		// if there's only 1 base, then it's pure.
		if (second==0) return 1;
		double result = (double) first / (double) (first + second);
		return result;
	}

	/**
	 * Only keep UMIs where the purity of the UMI is at least this threshold.
	 * @param threshold
	 */
	public void filterDataByUMIPurity (final double threshold) {
		List<String> removeKeys=new ArrayList<String>();

		for (String key: this.umiReadCounts.keySet()) {
			double purity = getUMIPurity(key);
			if (purity<threshold)
				removeKeys.add(key);
		}

		for (String key: removeKeys)
			this.umiReadCounts.remove(key);

	}

	public double getMeanUMIPurity () {
		List<String> umis = new ArrayList<String>(this.umis());

		double [] purity = new double [umis.size()];
		for (int i=0; i<purity.length; i++) {
			double d = getUMIPurity(umis.get(i));
			purity[i]=d;
		}

		double result = StatUtils.mean(purity);
		return result;
	}


	/**
	 * Get the most common and second most common alleles as the first two positions in the int array.
	 * If the reference base has been set for this object, then the this returns the count of the reference base, followed by the most common base that isn't the reference base.
	 * @param counts
	 * @return
	 */
	int [] getMostCommonCounts (final ObjectCounter<Character> counts, final Character referenceBase) {
		List<Character> keys = counts.getKeysOrderedByCount(true);

		int [] result=new int[2];
		Arrays.fill(result, 0);

		// the default behavior.
		if (referenceBase==null) {
			if (keys.size()>=1)
				result[0]=counts.getCountForKey(keys.get(0));
			if (keys.size()>=2)
				result[1]=counts.getCountForKey(keys.get(1));
		} else {
			result[0]=counts.getCountForKey(this.referenceBase);
			keys.remove(referenceBase);
			if (keys.size()>=1)
				result[1]=counts.getCountForKey(keys.get(0));

		}

		return (result);
	}



	/**
	 * Sets the reference base for this DAC.
	 * When this is set, BinomialStatistics use this base as the most common one
	 */
	public void setCommonBase (final char base) {
		this.referenceBase=new Character(base);
	}

	public char getReferenceBase () {
		return this.referenceBase;
	}



	/**
	 * Convert the pileup of bases and qualities into a filtered count of bases.
	 * @param p The input SNPUMIBasePileup to extract bases from.
	 * @return A count of bases that pass the baseQualityThreshold
	 */
	private ObjectCounter<Character> convertToObjectCounter(final SNPUMIBasePileup p) {
		ObjectCounter<Character> result = new ObjectCounter<Character>();
		// iterate over the base and quality.  They are the same size.
		List<Byte> quals = p.getQualities();
		Iterator<Character> bases = p.getBasesAsCharacters().iterator();
		for (Byte qual: quals) {
			Character base = bases.next();
			if (qual>=this.baseQualityThreshold) result.increment(base);
		}
		return result;
	}

	private void validateIncomingBasePileup (final SNPUMIBasePileup p) {
		// validate the gene/cell/snp/umi can be added.
		if (!this.gene.equals(p.getGene()) ||  !this.cell.equals(p.getCell()) || !this.snpInterval.equals(p.getSNPInterval()))
			throw new IllegalArgumentException("In DAC snp [" + this.snpInterval.getName() +"] gene [" + this.gene + "] cell [" + this.cell + "] Tried to add " + p.toString());
		if (this.umiReadCounts.containsKey(p.getMolecularBarcode()))
			throw new IllegalArgumentException("Already have UMI [" + p.getMolecularBarcode() +"] added to this DAC!");
	}

	/**
	 * Return a copy of this object, with the UMIs edit distance collapsed by Hamming distance.
	 * An edit distance of 0 is equivalent to a no-op.
	 * @param editDistance The threshold edit distance where UMIs are collapsed by Hamming distance.  UMIs separated by >= editDistance are collapsed.
	 * @return A shallow copy of the input DigitalAlleleCounts object with UMIs collapsed.
	 */
	public DigitalAlleleCounts collapseUMIs (final int editDistance) {
		DigitalAlleleCounts result = new DigitalAlleleCounts(this.snpInterval, this.gene, this.cell, this.baseQualityThreshold);
		// short circuit for ED=0.
		if (editDistance==0) {
			result.umiReadCounts = this.umiReadCounts;
			return (result);
		}
		// get the umis
		ObjectCounter<String> umis = getUMIsInReadOrder();
		// collapse all UMIs against each other.  "Primary" UMIs are keys.
		Map<String, List<String>> mapping = med.collapseBarcodes(umis, false, editDistance);

		// merge UMI results
		Map<String, ObjectCounter<Character>> map = new HashMap<String, ObjectCounter<Character>>();
		for (String umi: mapping.keySet()) {
			ObjectCounter<Character> r = new ObjectCounter<Character>(this.umiReadCounts.get(umi));
			for (String seconary: mapping.get(umi)) {
				ObjectCounter<Character> s = this.umiReadCounts.get(seconary);
				r.increment(s);
			}
			map.put(umi, r);
		}
		result.umiReadCounts=map;
		return (result);
	}

	/**
	 * Get a counter with the total number of reads for each UMI.
	 * @return
	 */
	private ObjectCounter<String> getUMIsInReadOrder () {
		ObjectCounter<String> result = new ObjectCounter<String>();
		for (String umi : this.umiReadCounts.keySet()) {
			int totalCount = this.umiReadCounts.get(umi).getTotalCount();
			result.incrementByCount(umi, totalCount);
		}
		return result;
	}

	public Interval getSnpInterval() {
		return snpInterval;
	}

	public String getGene() {
		return gene;
	}

	public String getCell() {
		return cell;
	}

	public int getBaseQualityThreshold() {
		return baseQualityThreshold;
	}

	/**
	 * Returns true if the object has no reads.
	 * @return
	 */
	public boolean isEmpty () {
		return this.getReadCounts().getTotalCount()==0;
	}

	@Override
	public String toString () {
		StringBuilder b = new StringBuilder();
		b.append("SNP [" + this.getSnpInterval()+"] ");
		b.append("Gene [" + this.getGene() +"] ");
		b.append("Cell [" + this.getCell() +"] ");
		b.append("Read counts [" + this.getReadCounts() +"] ");
		b.append("UMI counts [" + this.getUMIAlleleCount() +"] ");

		// only test if there are reads.
		if (!isEmpty()) {
			BinomialStatistics biReads = getReadConfidenceInterval(0.95);
			b.append("read pvalue [" + biReads.getBinomialPvalue() +"] ");
			b.append("read CI low [" + biReads.getBinomialConfidenceInterval().getLowerBound() + "] ");
			b.append("read CI high [" + biReads.getBinomialConfidenceInterval().getUpperBound() + "] ");
		}
		// only test if there are UMIs - this might collapse into if there are reads, then there are UMIs?
		// TODO should this fold into the same block as the reads binomial statistics?
		if (this.getUMIAlleleCount().getTotalCount()>0) {
			BinomialStatistics biUMIs = getUMIConfidenceInterval(0.95);
			b.append("UMI pvalue [" + biUMIs.getBinomialPvalue() +"] ");
			b.append("UMI CI low [" + biUMIs.getBinomialConfidenceInterval().getLowerBound() + "] ");
			b.append("UMI CI high [" + biUMIs.getBinomialConfidenceInterval().getUpperBound() + "] ");
		}
		return (b.toString());

	}


}
