package org.broadinstitute.dropseqrna.annotation;

import htsjdk.samtools.util.Interval;

public class GTFRecord implements Comparable<GTFRecord> {

	private Interval interval;
	private String geneID;
	private String geneName;
	private String transcriptName;
	private String transcriptID;
	private String transcriptType;
	private String featureType;
	
	
	//{chr, startPos, endPos, strand,geneName,geneID,transcriptName, transcriptID,transcriptType,featureType};
	public GTFRecord(String chromsome, int start, int end, boolean negativeStrand, String geneID, String geneName, String transcriptName, String transcriptID, String transcriptType, String featureType) {
		this.interval=new Interval(chromsome, start, end, negativeStrand, null);
		this.geneID=geneID;
		this.geneName=geneName;
		this.transcriptName=transcriptName;
		this.transcriptID=transcriptID;
		this.transcriptType=transcriptType;
		this.featureType=featureType;
	}
	
	public Interval getInterval () {
		return this.interval;
	}
	
	public String getStrandAsString() {
		if (interval.isNegativeStrand()) return ("-");
		return ("+");
	}
	
	public String getChromosome() {
		return this.interval.getSequence();
	}
	
	public int getStart() {
		return this.interval.getStart();
	}
	
	public int getEnd() {
		return this.interval.getEnd();
	}
	
	public boolean isNegativeStrand() {
		return this.interval.isNegativeStrand();
	}
	
	public String getGeneID() {
		return geneID;
	}

	public String getGeneName() {
		return geneName;
	}

	public String getTranscriptName() {
		return transcriptName;
	}

	public String getTranscriptID() {
		return transcriptID;
	}

	public String getTranscriptType() {
		return transcriptType;
	}

	public String getFeatureType() {
		return featureType;
	}
	
	@Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final GTFRecord that = (GTFRecord) o;
        if (!this.interval.equals(that.interval)) return false;
        if (!this.geneID.equals(that.geneID)) return false;
        if (!this.geneName.equals(that.geneName)) return false;
        if (!this.transcriptName.equals(that.transcriptName)) return false;
        if (!this.transcriptID.equals(that.transcriptID)) return false;
        if (this.transcriptType == null) {
            if (that.transcriptType != null) return false;
        } else if (!this.transcriptType.equals(that.transcriptType)) return false;
        if (!this.featureType.equals(that.featureType)) return false;
        
        return true;
    }

    @Override
    public int hashCode() {
        int result = interval.hashCode();
        result = 31 * result + geneID.hashCode();
        result = 31 * result + geneName.hashCode();
        result = 31 * result + transcriptName.hashCode();
        result = 31 * result + transcriptID.hashCode();
        if (transcriptType != null) {
            // E.g. ERCC
            result = 31 * result + transcriptType.hashCode();
        }
        result = 31 * result + featureType.hashCode();
        return result;
    }
    
    
   
    
   
	@Override
	public int compareTo(GTFRecord o) {
		return this.interval.compareTo(o.getInterval());
	}
	
	public String toString () {
		return (this.interval.toString() +" [" + this.geneName + " " + this.featureType+ "]");
	}
	
}
