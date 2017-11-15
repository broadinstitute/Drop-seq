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
