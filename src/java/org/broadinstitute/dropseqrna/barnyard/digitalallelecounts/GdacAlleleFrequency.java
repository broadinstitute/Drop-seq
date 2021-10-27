package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import org.apache.commons.lang.builder.CompareToBuilder;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

import htsjdk.samtools.util.Interval;

/**
 * Structure to store the allele frequency of a single variant.
 * @author nemesh
 *
 */
public class GdacAlleleFrequency implements Comparable<GdacAlleleFrequency> {

	private final Interval snpInterval;
	
	private final int refReadCount;
	private final int altReadCount;
	private final int refUmiCount;
	private final int altUmiCount;
	
	private final char refAllele;
	private final char altAllele;
	
	public GdacAlleleFrequency (DigitalAlleleCounts dac) {
		this.snpInterval=dac.getSnpInterval();
		this.refAllele=dac.getReferenceAllele();
		this.altAllele=dac.getAltAllele();
		
		ObjectCounter<Character> readCounts = dac.getReadCounts();
		ObjectCounter<Character> umiCounts = dac.getUMIAlleleCount();
		
		this.refReadCount=readCounts.getCountForKey(dac.getReferenceAllele());
		this.altReadCount=readCounts.getCountForKey(dac.getAltAllele());		
		this.refUmiCount=umiCounts.getCountForKey(dac.getReferenceAllele());
		this.altUmiCount=umiCounts.getCountForKey(dac.getAltAllele());
	}
	
	public GdacAlleleFrequency (Interval snpInterval, char refAllele, char altAllele, int refReadCount, int altReadCount, int refUmiCount, int altUmiCount) {
		this.snpInterval=snpInterval;
		this.refAllele=refAllele;
		this.altAllele=altAllele;
		this.refReadCount=refReadCount;
		this.altReadCount=altReadCount;
		this.refUmiCount=refUmiCount;
		this.altUmiCount=altUmiCount;
	}
	
	/**
	 * Merges these results with another GdacAlleleFrequency object.
	 * The interval, ref allele, and alt allele must match for both objcets.
	 * This sums the various counts and returns a new object with the results.
	 * @param other Another set of results to merge.  
	 * @return A new GdacAlleleFrequency object with merged counts from the current and other object.
	 */
	public GdacAlleleFrequency merge (GdacAlleleFrequency other) {
		if (this.refAllele!=other.getRefAllele())
			throw new IllegalArgumentException("Attempting to merge two GdacAlleleFrequency objects with different reference alleles");
		
		if (this.altAllele!=other.getAltAllele())
			throw new IllegalArgumentException("Attempting to merge two GdacAlleleFrequency objects with different alternate alleles");
		
		if (!this.snpInterval.equals(other.snpInterval))
			throw new IllegalArgumentException("Attempting to merge two GdacAlleleFrequency objects with different intervals");
		
		GdacAlleleFrequency result = new GdacAlleleFrequency(this.snpInterval, this.refAllele, this.altAllele, this.getRefReadCount()+other.getRefReadCount(),
				this.getAltReadCount()+other.getAltReadCount(), this.getRefUmiCount()+other.getRefUmiCount(), this.getAltUmiCount()+other.getAltUmiCount());
		
		return (result);				
	}
	
	public Interval getSnpInterval() {
		return snpInterval;
	}

	public int getRefReadCount() {
		return refReadCount;
	}

	public int getAltReadCount() {
		return altReadCount;
	}

	public int getRefUmiCount() {
		return refUmiCount;
	}

	public int getAltUmiCount() {
		return altUmiCount;
	}

	public char getRefAllele() {
		return refAllele;
	}

	public char getAltAllele() {
		return altAllele;
	}
	
	public Double getReadRatio () {
		if (this.refReadCount+this.altReadCount==0)
			return (null);
		double maf = (double) this.altReadCount/ (this.refReadCount+this.altReadCount);
		return (maf);
	}
	
	public Double getUMIRatio () {
		if (this.refUmiCount+this.altUmiCount==0)
			return (null);
		double maf = (double) this.altUmiCount/ (this.refUmiCount+this.altUmiCount);
		return (maf);
	}
	
	/**
	 * Get a unique string representation of the position in the genome.  Two GdacAlleleFrequency objects that have the same key are comparable.
	 * @return a unique string representation of the position in the genome
	 */
	public String getCompoundKey () {
		return this.snpInterval.toString()+":"+this.refAllele+":"+this.altAllele;		
	}
	
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + altAllele;
		result = prime * result + altReadCount;
		result = prime * result + altUmiCount;
		result = prime * result + refAllele;
		result = prime * result + refReadCount;
		result = prime * result + refUmiCount;
		result = prime * result + ((snpInterval == null) ? 0 : snpInterval.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		GdacAlleleFrequency other = (GdacAlleleFrequency) obj;
		if (altAllele != other.altAllele)
			return false;
		if (altReadCount != other.altReadCount)
			return false;
		if (altUmiCount != other.altUmiCount)
			return false;
		if (refAllele != other.refAllele)
			return false;
		if (refReadCount != other.refReadCount)
			return false;
		if (refUmiCount != other.refUmiCount)
			return false;
		if (snpInterval == null) {
			if (other.snpInterval != null)
				return false;
		} else if (!snpInterval.equals(other.snpInterval))
			return false;
		return true;
	}
	
	@Override
	public int compareTo(GdacAlleleFrequency o) {
		return new CompareToBuilder()	
					.append(this.snpInterval, o.snpInterval)
			       .append(this.refAllele, o.refAllele)
			       .append(this.altAllele, o.getAltAllele()).toComparison();				       			
	}
	
	

}
