/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.dropseqrna.annotation;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang.StringUtils;

/**
 * Holds one parsed line of a refFlat record.  Also good for writing them back to a file.
 * http://genome.ucsc.edu/goldenpath/gbdDescriptionsOld.html#RefFlat
 * 
 * Note that refFlat records are in 0-based half-open notation.  Data is converted from the assumed 1-based inclusive.  
 * That basically means you subtract one from the starts of exons/transcripts/CDS.
 * @author nemesh
 *
 */
public class RefFlatRecord {
	private String geneName;
	
	private String transcriptName;
	private String chromosome;
	private boolean isNegativeStrand;
	private int transcriptStart;
	private int transcriptEnd;
	private int cdsStart;
	private int cdsEnd;
	List<Integer> exonStarts;
	List<Integer> exonEnds;
	
	
	public RefFlatRecord (String geneName, String transcriptName, String chromosome, boolean isNegativeStrand, int transcriptStart, int transcriptEnd, 
			int cdsStart, int cdsEnd) {
		this.geneName=geneName;
		this.transcriptName=transcriptName;
		this.chromosome=chromosome;
		this.isNegativeStrand=isNegativeStrand;
		this.transcriptStart=transcriptStart-1;
		this.transcriptEnd=transcriptEnd;
		this.cdsStart=cdsStart-1;
		this.cdsEnd=cdsEnd;
		this.exonStarts = new ArrayList<Integer>();
		this.exonEnds = new ArrayList<Integer>();	
	}
	
	public String getGeneName() {
		return geneName;
	}

	public void setGeneName(String geneName) {
		this.geneName = geneName;
	}

	public String getTranscriptName() {
		return transcriptName;
	}

	public String getChromosome() {
		return chromosome;
	}

	public boolean isNegativeStrand() {
		return isNegativeStrand;
	}

	public int getTranscriptStart() {
		return transcriptStart;
	}

	public int getTranscriptEnd() {
		return transcriptEnd;
	}

	public int getCdsStart() {
		return cdsStart;
	}

	public int getCdsEnd() {
		return cdsEnd;
	}

	public List<Integer> getExonStarts() {
		return exonStarts;
	}
	
	public List<Integer> getExonEnds() {
		return exonEnds;
	}
	
	public void addExonStart (int start) {
		this.exonStarts.add(start-1);
	}
	
	public void addExonEnd (int end) {
		this.exonEnds.add(end);
	}
	
	public String getFormattedExonString(List<Integer> positions) {
		StringBuilder b = new StringBuilder();
		for (Integer i: positions) {
			b.append(i.toString());
			b.append(",");
		}
		return (b.toString());
	}
	
	/**
	 * This is the tab formatted output expected for a RefFlat file.
	 */
	public String toString() {
		List<String> r = new ArrayList<String>(11);
		r.add(this.geneName);
		r.add(this.transcriptName);
		r.add(this.chromosome);
		r.add(AnnotationUtils.strandToString(!this.isNegativeStrand));
		r.add(Integer.toString(this.transcriptStart));
		r.add(Integer.toString(this.transcriptEnd));
		r.add(Integer.toString(this.cdsStart));
		r.add(Integer.toString(this.cdsEnd));
		r.add(Integer.toString(this.getExonStarts().size()));
		r.add(getFormattedExonString(this.exonStarts));
		r.add(getFormattedExonString(this.exonEnds));
		String h = StringUtils.join(r, "\t");
		return (h);
	}
	
	@Override
    public boolean equals(final Object o) {
		 if (this == o) return true;
         if (o == null || getClass() != o.getClass()) return false;

         final RefFlatRecord that = (RefFlatRecord) o;
         if (!this.geneName.equals(that.geneName)) return false;
         if (!this.transcriptName.equals(that.transcriptName)) return false;
         if (!this.chromosome.equals(that.chromosome)) return false;
         if (this.transcriptStart!=that.transcriptStart) return false;
         if (this.transcriptEnd!=that.transcriptEnd) return false;
         if (this.cdsStart!=that.cdsStart) return false;
         if (this.cdsEnd!=that.cdsEnd) return false;
         if (!getFormattedExonString(this.exonStarts).equals(getFormattedExonString(that.exonStarts))) return false;
         if (!getFormattedExonString(this.exonEnds).equals(getFormattedExonString(that.exonEnds))) return false;
         return true;
	}
	
}
