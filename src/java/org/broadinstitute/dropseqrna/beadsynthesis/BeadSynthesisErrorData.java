package org.broadinstitute.dropseqrna.beadsynthesis;

import java.util.Collection;
import java.util.List;

import org.apache.commons.lang.math.NumberUtils;
import org.broadinstitute.dropseqrna.utils.BaseDistributionMetric;
import org.broadinstitute.dropseqrna.utils.BaseDistributionMetricCollection;
import org.broadinstitute.dropseqrna.utils.Bases;

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
	 * Classify the error type if any for this data.
	 * @param threshold what fraction of the bases at a position must be all one base.  
	 * @return The error type for the data.
	 */
	public BeadSynthesisErrorTypes getErrorType (double threshold) {
		int errorPosition = getErrorBase(threshold);
		int polyTPos = getPolyTErrorPosition(threshold);		
		if (errorPosition==polyTPos & errorPosition!=-1) {
			return BeadSynthesisErrorTypes.POLY_T_ERROR;				
		}
		if (hasSingleUMIError(threshold)) {
			return BeadSynthesisErrorTypes.SINGLE_UMI;
		}
		if (errorPosition!=polyTPos) {
			return BeadSynthesisErrorTypes.OTHER_ERROR;
		}
		return BeadSynthesisErrorTypes.NO_ERROR;
	}
	
	/**
	 * A special case error, the polyT error occurs when the dominant base
	 * is T across 1 or more positions. 
	 * @param threshold what fraction of the bases at this position must be T.
	 * @return The base in the UMI where the error begins.  If this is the last base of an length 8 umi, the result would be 8.
	 * If no base position is predominantly polyT, return -1;
	 * IE: the return is 1 based.
	 */
	public int getPolyTErrorPosition (double threshold) {
		double [] freq = getPolyTFrequency();
		for (int position=0; position<freq.length; position++) { 
			if (freq[position]>=threshold) return (position+1);
		}
		return -1;
	}
	
	
	
	/**
	 * A special case error where all bases appear fixed.
	 * @param threshold what fraction of a base /sum(bases) must be fixed at all positions
	 * @return if all bases are skewed to represent <threshold> fraction of bases at all positions return true.  Otherwise, return false.
	 */
	public boolean hasSingleUMIError (double threshold) {
		double [] data = synthesisErrorMetric();
		for (double d: data) {
			if (d<threshold) return (false);
		}
		return (true);
	}
	
	/**
	 * A special case error, the polyT error occurs when the dominant base
	 * Get the frequency of T at each position
	 */
	private double [] getPolyTFrequency () {
		char base = Bases.T.getBase();
		List<Integer> basePositions = baseCounts.getPositions();
		double [] result = new double [basePositions.size()];
		for (int position=0; position<result.length; position++) {
			BaseDistributionMetric bdm =this.baseCounts.getDistributionAtPosition(position);
			double freq = (double) bdm.getCount(base) / (double) bdm.getTotalCount();
			result[position]=freq;
		}
		return (result);
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
