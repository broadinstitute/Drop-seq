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

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.samtools.util.ProgressLogger;
import picard.annotation.AnnotationException;
import picard.annotation.Gene.Transcript.Exon;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class EnhanceGTFRecords {

    private static final Log LOG = Log.getInstance(EnhanceGTFRecords.class);
    public static final String EXON_FEATURE_TYPE = "exon";
    public static final String INTRON_FEATURE_TYPE = "intron";
    public static final String CONSENSUS_INTRON_FEATURE_TYPE = "consensus_intron";
    public static final String GENE_FEATURE_TYPE = "gene";
    public static final String TRANSCRIPT_FEATURE_TYPE = "transcript";


    /**
	 * @param records GTFRecords, potentially from many genes
	 * @return The input records, plus gene, transcript, intron and conserved_intron records.
	 */
	public List<GTFRecord> enhanceGTFRecords (Iterator<GTFRecord> records) {
        final ProgressLogger progressLogger = new ProgressLogger(LOG, 100, "enhancing", "genes");
        List<GTFRecord> result = new ArrayList<>();
        final GeneFromGTFBuilder geneBuilder = new GeneFromGTFBuilder(records);
        while (geneBuilder.hasNext()) {
            try {
                progressLogger.record();
                GeneFromGTF gene = geneBuilder.next();
                progressLogger.record(gene.getName(), 0);
                LOG.debug("Enhancing Gene [" + gene.getName() + "]");
                result.addAll(enhanceGene(gene));
            } catch (AnnotationException e) {
                LOG.info(e.getMessage() + " -- skipping");
            }
        }
		return result;
	}
	
	public List<GTFRecord> enhanceGene (GeneFromGTF g) {
		List<GTFRecord> result = new ArrayList<>();
		// add gene entry
		result.add(getGTFRecord(g));
		//add transcript entries.
		for (GeneFromGTF.TranscriptFromGTF t: g.getTranscripts()) {
			List<GTFRecord> records = enhanceTranscript(t);
			result.addAll(records);
		}
        addConsensusIntrons(result);
		return (result);
	}

    private void addConsensusIntrons(final List<GTFRecord> records) {
        final OverlapDetector<GTFRecord> exons = new OverlapDetector<>(0, 0);
        for (final GTFRecord record: records) {
            if (record.getFeatureType().equals(EXON_FEATURE_TYPE)) {
                exons.addLhs(record, record.getInterval());
            }
        }
        final ArrayList<GTFRecord> consensusIntrons = new ArrayList<>();
        for (final GTFRecord record: records) {
            if (record.getFeatureType().equals(INTRON_FEATURE_TYPE) && exons.getOverlaps(record.getInterval()).isEmpty()) {
                consensusIntrons.add(new GTFRecord(
                        record.getChromosome(),
                        record.getStart(),
                        record.getEnd(),
                        record.isNegativeStrand(),
                        record.getGeneID(),
                        record.getGeneName(),
                        record.getTranscriptName(),
                        record.getTranscriptID(),
                        record.getTranscriptType(),
                        CONSENSUS_INTRON_FEATURE_TYPE,
                        record.getGeneVersion()));
            }
        }
        records.addAll(consensusIntrons);
    }

	private GTFRecord getGTFRecord(GeneFromGTF g) {
        //noinspection deprecation
        return new GTFRecord(g.getContig(), g.getStart(), g.getEnd(), g.isNegativeStrand(),
				g.getGeneID(), g.getName(), null, null, g.getTranscriptType(), GENE_FEATURE_TYPE, g.getGeneVersion());
	}
	
	
	private List<GTFRecord> enhanceTranscript (GeneFromGTF.TranscriptFromGTF t) {
		List<GTFRecord> result = new ArrayList<>();
		GeneFromGTF g = t.getGene();
		@SuppressWarnings("deprecation") GTFRecord transcriptRecord = new GTFRecord(g.getContig(), t.codingStart, t.codingEnd, g.isNegativeStrand(),
                g.getGeneID(), g.getName(), t.getTranscriptName(), t.getTranscriptID(), t.getTranscriptType(),
                TRANSCRIPT_FEATURE_TYPE, g.getGeneVersion());
		result.add(transcriptRecord);
		
		List<Interval> introns = getIntronIntervals(getIntervals(t.exons));
		List<GTFRecord> intronR = getGTFRecordsFromIntronIntervals(introns, t);
		result.addAll(intronR);
		List<GTFRecord> exons = getGTFRecordsFromExons(t);
		result.addAll(exons);
		return (result);
		
	}
	
	
	private List<GTFRecord> getGTFRecordsFromExons(GeneFromGTF.TranscriptFromGTF t) {
		List<GTFRecord> result = new ArrayList<>();
		GeneFromGTF g = t.getGene();
		for (Exon e: t.exons) {
			@SuppressWarnings("deprecation") GTFRecord exonRecord = new GTFRecord(g.getContig(), e.start, e.end, g.isNegativeStrand(), g.getGeneID(),
                    g.getName(), t.getTranscriptName(), t.getTranscriptID(), t.getTranscriptType(), EXON_FEATURE_TYPE,
                    g.getGeneVersion());
			result.add(exonRecord);
		}
		return (result);
	}
	
	private List<GTFRecord> getGTFRecordsFromIntronIntervals (List<Interval> introns, GeneFromGTF.TranscriptFromGTF t) {
		List<GTFRecord> result = new ArrayList<>();
		GeneFromGTF g = t.getGene();
		for (Interval i: introns) {
			@SuppressWarnings("deprecation") GTFRecord intronRecord = new GTFRecord(g.getContig(), i.getStart(), i.getEnd(), g.isNegativeStrand(),
                    g.getGeneID(), g.getName(), t.getTranscriptName(), t.getTranscriptID(), t.getTranscriptType(),
                    INTRON_FEATURE_TYPE, g.getGeneVersion());
			result.add(intronRecord);
		}		
		return (result);
	}
	
	
	List<Interval> getIntronIntervals(List<Interval> exons) {
		if (exons.size()==1) return new ArrayList<>(0);
		
		List<Interval> result = new ArrayList<>(exons.size()-1);
		for (int i=0; i<exons.size()-1; i++) {
			int s = exons.get(i).getEnd()+1;
			int e = exons.get(i+1).getStart()-1;
			Interval intron = new Interval(null, s, e);
			result.add(intron);
		}
		return (result);
	}
	
	List<Interval> getIntervals (Exon [] exons) {
		List<Interval> result = new ArrayList<>();
		for (Exon e: exons) {
			result.add(new Interval(null, e.start, e.end));
		}
		return result;
	}
}
