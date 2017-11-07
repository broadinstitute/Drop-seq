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

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import picard.annotation.AnnotationException;
import picard.annotation.Gene;


@SuppressWarnings("Convert2Diamond")
public class GeneFromGTF extends Gene implements Iterable<Gene.Transcript> {
	
	private final String geneID;
	private final String transcriptType;
	private final String featureType;
    private final Integer geneVersion;
	
	private final Map<String, TranscriptFromGTF> transcripts = new HashMap<String, TranscriptFromGTF>();
	
	public GeneFromGTF(String sequence, int start, int end, boolean negative, String name, String featureType,
                       String geneID, String transcriptType, Integer geneVersion) {
		super(sequence, start, end, negative, name);
		this.featureType=featureType;
		this.geneID=geneID;
		this.transcriptType=transcriptType;
        this.geneVersion = geneVersion;
	}
	
	
	public String getGeneID() {
		return geneID;
	}

	public String getTranscriptType() {
		return transcriptType;
	}

	public String getFeatureType() {
		return featureType;
	}

    public Integer getGeneVersion() { return geneVersion; }

    @SuppressWarnings({ "rawtypes", "unchecked" })
	@Override
	public Iterator<Transcript> iterator() {
        return (Iterator) transcripts.values().iterator();
    }
	
	public Collection<TranscriptFromGTF> getTranscripts() {
		return transcripts.values();
	}
	
	public TranscriptFromGTF addTranscript(final String name, final int transcriptionStart, final int transcriptionEnd, final int codingStart, final int codingEnd, 
			final int numExons, String transcriptName, String transcriptID, String transcriptType) {
        if (transcripts.containsKey(name)) {
            throw new AnnotationException("Transcript " + name + " for gene " + this.getName() + " appears more than once");
        }
        else {
            final TranscriptFromGTF tx = new TranscriptFromGTF(name, transcriptionStart, transcriptionEnd, codingStart, codingEnd, numExons, 
            		transcriptName, transcriptID,  transcriptType);
            transcripts.put(name, tx);
            return tx;
        }
    }
	
	 @Override
     public boolean equals(final Object o) {
         if (this == o) return true;
         if (o == null || getClass() != o.getClass()) return false;

         final GeneFromGTF that = (GeneFromGTF) o;

         if (this.getStart() != that.getStart()) return false;
         if (this.getEnd() != that.getEnd()) return false;
         if (!this.getName().equals(that.getName())) return false;
         if (!this.getContig().equals(that.getContig())) return false;
         if (!this.getGeneID().equals(that.getGeneID())) return false;

         return true;
     }

     @Override
     public int hashCode() {
         int result = this.getContig().hashCode();
         result = 31 * result + this.getStart();
         result = 31 * result + this.getEnd();
         result = 31 * result + this.getName().hashCode();
         return result;
     }
			
	public class TranscriptFromGTF extends Gene.Transcript {

		private final String transcriptName;
		private final String transcriptID;
		private final String transcriptType;
		
		public TranscriptFromGTF(String name, int transcriptionStart, int transcriptionEnd, int codingStart, int codingEnd,
				int numExons, String transcriptName, String transcriptID, String transcriptType) {
			super(name, transcriptionStart, transcriptionEnd, codingStart, codingEnd,
					numExons);
			this.transcriptName=transcriptName;
			this.transcriptID=transcriptID;
			this.transcriptType=transcriptType;
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
		
		public GeneFromGTF getGene() {
	        return GeneFromGTF.this;
	    }
		
		
		
	}

	
	
	
}
