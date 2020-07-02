package org.broadinstitute.dropseqrna.censusseq;

import org.apache.commons.math3.stat.descriptive.moment.Mean;

import htsjdk.samtools.metrics.MetricBase;

public class CsiMetrics extends MetricBase {
	public int NUM_SNPS;
	public int REF_COUNT;
	public int ALT_COUNT;
	public double REF_DOSAGE;
	public double ALT_DOSAGE;
	public double MEAN_MAF;
	public double FRAC_ALT;
	public double FRAC_DOSAGE_ALT;
	public double SEQUENCING_ERROR_RATE = 0.001;
	public double PCT_CONTAMINATION;
	public double OBSERVED_SEQUENCING_ERROR_RATE;
	// public double PCT_CONTAMINATION_EMPIRIC_ERROR_RATE;
	// public double PCT_CONTAMINATION_BQ_AWARE;
	public String KNOWN_CONTAMINATION_RATE;
	
	private Mean afMean = new Mean();
	private Mean errorProbMean = new Mean();

	public void addAlleleFreq(final double alleleFreq) {
		this.afMean.increment(alleleFreq);
	}

	public void addBaseErrorProbability(final double errorRate) {
		errorProbMean.increment(errorRate);
	}
	
	public double getAltFrequencyByCount () {		
		return (double) this.ALT_COUNT / (double) (this.ALT_COUNT + this.REF_COUNT);
	}
	
	public double getAltFreqByDosage () {
		return (double) this.ALT_DOSAGE / (double) (this.ALT_DOSAGE + this.REF_DOSAGE);
	}
	
	public double getObservedSequencingErrorRate () {
		return errorProbMean.getResult();
	}

	public void calculateStats() {
		MEAN_MAF = this.afMean.getResult();
		OBSERVED_SEQUENCING_ERROR_RATE = getObservedSequencingErrorRate();
		this.FRAC_ALT = getAltFrequencyByCount();			
		this.PCT_CONTAMINATION = ((this.FRAC_ALT - this.SEQUENCING_ERROR_RATE) / this.MEAN_MAF) * 100;
		//this.PCT_CONTAMINATION_EMPIRIC_ERROR_RATE = ((this.FRAC_ALT - this.OBSERVED_SEQUENCING_ERROR_RATE)
		// 		/ this.MEAN_MAF) * 100;
		// if (this.PCT_CONTAMINATION_EMPIRIC_ERROR_RATE<0) this.PCT_CONTAMINATION_EMPIRIC_ERROR_RATE=0;
		if (this.PCT_CONTAMINATION<0) this.PCT_CONTAMINATION=0;
		
		// base quality aware [experimental]
		this.FRAC_DOSAGE_ALT= getAltFreqByDosage();			
		// this.PCT_CONTAMINATION_BQ_AWARE = ((this.FRAC_DOSAGE_ALT - this.OBSERVED_SEQUENCING_ERROR_RATE)
		// 		/ this.MEAN_MAF) * 100;
		
	}
}
