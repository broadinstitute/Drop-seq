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
package org.broadinstitute.dropseqrna.beadsynthesis;

import java.util.Collection;
import java.util.Comparator;
import java.util.List;

import org.apache.commons.lang.math.NumberUtils;
import org.broadinstitute.dropseqrna.utils.BaseDistributionMetric;
import org.broadinstitute.dropseqrna.utils.BaseDistributionMetricCollection;
import org.broadinstitute.dropseqrna.utils.Bases;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

public class BeadSynthesisErrorData {

	private final String cellBarcode;
	private BaseDistributionMetricCollection baseCounts;
	private ObjectCounter<String> umiCounts;
	private int numReads;
	private int numTranscripts;

	// cached results
	private boolean dataChanged;
	private BeadSynthesisErrorTypes errorTypeCached=null;
	private BeadSynthesisErrorTypes errorTypeExtendedCached=null;

	// the number of unique sequences without collapse.
	private int numUMIs;

	public BeadSynthesisErrorData (final String cellBarcode) {
		this.cellBarcode=cellBarcode;
		this.baseCounts = new BaseDistributionMetricCollection();
		this.umiCounts = new ObjectCounter<>();
		this.dataChanged=true;
		this.numReads=0;
		this.numTranscripts=0;
		this.numUMIs=-1;
	}

	/**
	 * Build an object stub for testing that has the cell barcode and the error type.
	 * This error type is set for all calls to getErrorType.
	 * @param cellBarcode
	 * @param errorType
	 * @return
	 */
	static BeadSynthesisErrorData getInstance(final String cellBarcode, final BeadSynthesisErrorTypes errorType, final BaseDistributionMetricCollection baseCounts) {
		BeadSynthesisErrorData r = new BeadSynthesisErrorData(cellBarcode);
		r.baseCounts=baseCounts;
		r.errorTypeCached=errorType;
		r.errorTypeExtendedCached=errorType;
		r.dataChanged=false;
		r.numReads=0;
		r.numTranscripts=0;
		r.numUMIs=0;
		return null;
	}

	//maybe add a finalize() method that gets the count of the total UMIs in the object counter, caches the number, and throws away the object counter since that might be expensive
	// to store.
	// if you finalize, you also need to run AbstractDetectBeadSynthesisErrors.getEnhancedErrorType first, so that should probably be part of this class, and not an "enhanced" error.

	public void addUMI (final String umi) {
		//umiCounts++;
		this.umiCounts.increment(umi);
		baseCounts.addBases(umi);
		this.dataChanged=true;
	}

	public void addUMI (final Collection <String> umis) {
		for (String umi: umis)
			addUMI(umi);
	}

	/**
	 * No more data can be added after this is called.  The collection of UMIs will be discarded to save space.
	 * Whatever call was made to getErrorType will be cached, so one of the methods must be called before finalize.
	 */
	@Override
	public void finalize () {
		if (this.errorTypeCached==null && this.errorTypeExtendedCached==null)
			throw new IllegalStateException("Can't finalize until call to getErrorType is made!");
		// cache the number of UMIs.
		getUMICount();
		// remove the ObjectCounter to save memory
		this.umiCounts=null;
	}

	public void incrementReads (final int numReads) {
		this.numReads+=numReads;
	}

	public void incrementTranscripts(final int numTranscripts) {
		this.numTranscripts+=numTranscripts;
	}

	public int getNumReads () {
		return this.numReads;
	}

	public int getNumTranscripts() {
		return this.numTranscripts;
	}

	public int getUMICount() {
		if (this.numUMIs==-1 || this.dataChanged)
			this.numUMIs=this.umiCounts.getTotalCount();
		return this.numUMIs;
	}

	public ObjectCounter<String> getUMICounts() {
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
	// the caching is ugly and needs to be cleaned up.
	public BeadSynthesisErrorTypes getErrorType (final double threshold) {
		// return cached result
		if (!this.dataChanged & this.errorTypeCached!=null) return (this.errorTypeCached);
		this.dataChanged=false;

		int errorPosition = getErrorBase(threshold);
		int polyTPos = getPolyTErrorPosition(threshold);
		if (errorPosition==polyTPos & errorPosition!=-1) {
			BeadSynthesisErrorTypes t = BeadSynthesisErrorTypes.SYNTH_MISSING_BASE;
			this.errorTypeCached=t;
			return t;
		}
		else if (hasSingleUMIError(threshold)) {
			BeadSynthesisErrorTypes t = BeadSynthesisErrorTypes.SINGLE_UMI;
			this.errorTypeCached=t;
			return t;
		}
		else if (hasFixedFirstBase(threshold)) {
			BeadSynthesisErrorTypes t = BeadSynthesisErrorTypes.FIXED_FIRST_BASE;
			this.errorTypeCached=t;
			return t;
		}
		else if (errorPosition!=polyTPos) {
			BeadSynthesisErrorTypes t = BeadSynthesisErrorTypes.OTHER_ERROR;
			this.errorTypeCached=t;
			return t;
		}
		BeadSynthesisErrorTypes t = BeadSynthesisErrorTypes.NO_ERROR;
		this.errorTypeCached=t;
		return t;
	}

	BeadSynthesisErrorTypes getErrorType (final double extremeBaseRatio, final DetectPrimerInUMI detectPrimerTool, final Integer editDistanceToPrimer) {
		if (!this.dataChanged & this.errorTypeExtendedCached!=null) return (this.errorTypeExtendedCached);
		this.dataChanged=false;

		BeadSynthesisErrorTypes errorType = this.getErrorType(extremeBaseRatio);
		// cache the extended type as the regular type, but there's a chance to change if the detectPrimerTool says so.
		this.errorTypeExtendedCached=this.errorTypeCached;

		//base case, error is not a single UMI.
		if (errorType!=BeadSynthesisErrorTypes.SINGLE_UMI)
			return this.errorTypeExtendedCached;
		// if there's a primer, run detection.
		if (detectPrimerTool!=null & editDistanceToPrimer!=null) {
			// a single UMI-style error, does the most common UMI match the primer?
			String singleUMI = this.getUMICounts().getKeysOrderedByCount(true).get(0);
			boolean primerDetected = detectPrimerTool.isStringInPrimer(singleUMI, editDistanceToPrimer);
			if (primerDetected)
				this.errorTypeExtendedCached=BeadSynthesisErrorTypes.PRIMER;
		}
		return this.errorTypeExtendedCached;
	}

	/**
	 * A special case error, the polyT error occurs when the dominant base
	 * is T across 1 or more positions. The first base observed must be an 8, and subsequent bases must be adjacent
	 * @param threshold what fraction of the bases at this position must be T.
	 * @return The base in the UMI where the error begins.  If this is the last base of an length 8 umi, the result would be 8.
	 * If no base position is predominantly polyT, return -1;
	 * IE: the return is 1 based.
	 */
	public int getPolyTErrorPosition (final double threshold) {
		double [] freq = getPolyTFrequency();

		int errorBase=-2;
		for (int position=(freq.length)-1; position>=0; position--) {
			double v = freq[position];
			if (v>=threshold)
				errorBase=position;
			else
				break;
		}
		return errorBase+1;
	}



	/**
	 * A special case error where all bases appear fixed.
	 * @param threshold what fraction of a base /sum(bases) must be fixed at all positions
	 * @return if all bases are skewed to represent <threshold> fraction of bases at all positions return true.  Otherwise, return false.
	 */
	public boolean hasSingleUMIError (final double threshold) {
		double [] data = synthesisErrorMetric();
		for (double d: data)
			if (d<threshold) return (false);
		return (true);
	}

	/**
	 * Only the first base of the UMI is fixed > threshold
	 * @param threshold
	 * @return
	 */
	public boolean hasFixedFirstBase (final double threshold) {
		double [] data = synthesisErrorMetric();
		boolean pos0Fixed=data[0]>= threshold;
		// all positions after the first position
		boolean otherPositionsNotFixed=true;
		for (int i=1; i<data.length; i++)
			if (data[i]>threshold) {
				otherPositionsNotFixed=false;
				break;
			}
		boolean fixedFirstBase=(pos0Fixed &  otherPositionsNotFixed);
		return fixedFirstBase;
	}

	/**
	 * A special case error, the polyT error occurs when the dominant base
	 * Get the frequency of T at each position
	 */
	public double [] getPolyTFrequency () {
		//if (!this.dataChanged & this.polyTFreq!=null) return (this.polyTFreq);
		this.dataChanged=false;

		char base = Bases.T.getBase();
		List<Integer> basePositions = baseCounts.getPositions();
		double [] result = new double [basePositions.size()];
		for (int position=0; position<result.length; position++) {
			BaseDistributionMetric bdm =this.baseCounts.getDistributionAtPosition(position);
			double freq = (double) bdm.getCount(base) / (double) bdm.getTotalCount();
			result[position]=freq;
		}
		// cache results if you needed to compute.
		//this.polyTFreq=result;

		return (result);
	}

	/**
	 * Is the cell polyT biased at a particular base location?
	 * @param position The position in the UMI (1 based)
	 * @param threshold If the UMI at least this fraction T, the cell barcode is biased.
	 * @return
	 */
	public boolean isPolyTBiasedAtPosition(final int position, final double threshold) {
		double [] freq = getPolyTFrequency();
		if (position<1 | position > freq.length+1)
			throw new IllegalArgumentException("Requested position ["+ position+"] is outside the expected range of the UMI available bases 1-"+ freq.length+1);
		if (freq[position-1]>=threshold) return true;
		return false;
	}



	/**
	 * Convenience metric - return if the pileup of UMIs has a synthesis error.
	 * This gets the error metric for each base, finds the highest value, and compares it to the threshold.
	 * If the result is higher than the threshold return true.
	 * @return If the result is higher than the threshold return true.
	 */
	public boolean hasSynthesisError (final double threshold) {
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
	public int getErrorBase (final double threshold) {
		double [] m = synthesisErrorMetric();
		for (int i=0; i<m.length; i++)
			if (m[i]>=threshold) return (i+1);
		return -1;
	}

	/**
	 * Check if any base position of the UMI is dominated by a single base-type (A/C/G/T).
	 * Return the frequency of the most commonly occurring base observed at each position along the UMI
	 * @return
	 */
	public double [] synthesisErrorMetric() {
		//if (!this.dataChanged & this.synthesisErrorMetric!=null) return (this.synthesisErrorMetric);
		this.dataChanged=false;
		List<Integer> basePositions = baseCounts.getPositions();
		double [] result = new double [basePositions.size()];

		for (int i: basePositions)
			result[i]=getMostCommonBaseFrequency(i);
		//this.synthesisErrorMetric=result;
		return result;
	}

	/**
	 *
	 * @param position
	 * @return
	 */
	private double getMostCommonBaseFrequency (final int position) {
		BaseDistributionMetric bdm =this.baseCounts.getDistributionAtPosition(position);
		// cast as double once to avoid doing it over and over.
		double totalCount = bdm.getTotalCount();
		double maxFreq=0;

		for (Bases b: Bases.values()) {
			char bb = b.getBase();
			int count = bdm.getCount(bb);
			double freq = count / totalCount;
			if (freq>maxFreq)
				maxFreq=freq;
		}
		return (maxFreq);
	}

	public static class SizeComparator implements Comparator<BeadSynthesisErrorData> {
        @Override
        public int compare(final BeadSynthesisErrorData d1, final BeadSynthesisErrorData d2) {
        	int cmp = Integer.compare(d2.getUMICount(), d1.getUMICount());
            if (cmp != 0)
    			return cmp;
            return d1.getCellBarcode().compareTo(d2.getCellBarcode());
        }
    }

}
