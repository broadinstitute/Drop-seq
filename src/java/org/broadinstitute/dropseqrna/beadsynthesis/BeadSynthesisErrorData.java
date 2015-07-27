package org.broadinstitute.dropseqrna.beadsynthesis;

import java.util.Collection;
import java.util.List;

import org.apache.commons.lang.math.NumberUtils;
import org.broadinstitute.dropseqrna.utils.BaseDistributionMetric;
import org.broadinstitute.dropseqrna.utils.BaseDistributionMetricCollection;
import org.broadinstitute.dropseqrna.utils.Bases;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

public class BeadSynthesisErrorData {

	private final String cellBarcode;
	private BaseDistributionMetricCollection baseCounts;
	// private ObjectCounter<String> umiCounts;
	private int umiCounts;
	
	public BeadSynthesisErrorData (String cellBarcode) {
		this.cellBarcode=cellBarcode;
		this.baseCounts = new BaseDistributionMetricCollection();
		//this.umiCounts = new ObjectCounter<String>();
		umiCounts=0;
	}
	
	public void addUMI (String umi) {
		umiCounts++;
		//this.umiCounts.increment(umi);
		baseCounts.addBases(umi);
	}
	
	public void addUMI (Collection <String> umis) {
		for (String umi: umis) {
			addUMI(umi);
		}
	}
	
	public int getUMICount() {
		// return this.umiCounts.getSize();
		return this.umiCounts;
	}
	
	public String getCellBarcode() {
		return this.cellBarcode;
	}
	
	public BaseDistributionMetricCollection getBaseCounts () {
		return this.baseCounts;
	}
	
	public int getBaseLength () {
		return this.baseCounts.getPositions().size();
	}
	
	
	
	/**
	 * Convenience metric - return if the pileup of UMIs has a synthesis error.
	 * This gets the error metric for each base, finds the highest value, and compares it to the threshold.
	 * If the result is higher than the threshold return true.
	 * @return If the result is higher than the threshold return true.
	 */
	public boolean hasSynthesisError (double threshold) {
		double [] m = synthesisErrorMetric();
		double max = NumberUtils.max(m);
		return (max >= threshold);
	}
	
	/**
	 * Get the first base in the UMI that exhibits a synthesis error - the fraction of the 
	 * most common base at this position is >= the threshold.
	 * @param threshold The fraction of the offending base compared to all base counts at the position.  From 0-1.
	 * @return The position of the error, or -1 for no error.  This return is 1 based.
	 */
	public int getErrorBase (double threshold) {
		double [] m = synthesisErrorMetric();
		for (int i=0; i<m.length; i++) {
			if (m[i]>=threshold) return (i+1);
		}
		return -1;
	}
	
	/**
	 * Check if any base position of the UMI is dominated by a single base-type (A/C/G/T).
	 * Return the frequency of the most commonly occurring base observed at each position along the UMI  
	 * @return
	 */
	public double [] synthesisErrorMetric() {
		List<Integer> basePositions = baseCounts.getPositions();
		double [] result = new double [basePositions.size()];
		
		for (int i: basePositions) {
			result[i]=getMostCommonBaseFrequency(i);
		}
		return result;
	}
	
	/**
	 * 
	 * @param position
	 * @return
	 */
	private double getMostCommonBaseFrequency (int position) {
		BaseDistributionMetric bdm =this.baseCounts.getDistributionAtPosition(position);
		// cast as double once to avoid doing it over and over.
		double totalCount = (double) bdm.getTotalCount();
		double maxFreq=0;
		
		for (Bases b: Bases.values()) {
			char bb = b.getBase();
			int count = bdm.getCount(bb);
			double freq = (double) count / totalCount;
			if (freq>maxFreq) {
				maxFreq=freq;
			}
		}
		
		return (maxFreq);
		 
	}
	
	
	
}
