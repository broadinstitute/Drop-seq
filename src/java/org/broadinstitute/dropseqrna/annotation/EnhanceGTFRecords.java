package org.broadinstitute.dropseqrna.annotation;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;
import picard.annotation.Gene.Transcript.Exon;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class EnhanceGTFRecords {

    private static final Log log = Log.getInstance(ReduceGTF.class);
    public static final String EXON_FEATURE_TYPE = "exon";
    public static final String INTRON_FEATURE_TYPE = "intron";
    public static final String CONSENSUS_INTRON_FEATURE_TYPE = "consensus_intron";
    public static final String GENE_FEATURE_TYPE = "gene";
    public static final String TRANSCRIPT_FEATURE_TYPE = "transcript";


    /**
	 * For a list of GTFRecord objects for a single gene, add the gene,transcript,intron, and conserved intron GTFRecords, which can be determined by the GTFRecords.
	 * @param records
	 * @return
	 */
	public List<GTFRecord> enhanceGTFRecords (List<GTFRecord> records) {
		OverlapDetector<GeneFromGTF> od =getOverlapDetector(records);
		List<GTFRecord> result = enhanceGTFRecords(od);
		return (result);
	}
	
	public List<GTFRecord> enhanceGTFRecords (OverlapDetector<GeneFromGTF> od) {
		List<GTFRecord> result = new ArrayList<GTFRecord> ();
		for (GeneFromGTF g: od.getAll()) {
			
			log.info("Enhancing Gene [" + g.getName() +"]");
			List<GTFRecord> t = enhanceGene(g);
			result.addAll(t);
		}
		return (result);
	}
	
	public List<GTFRecord> enhanceGene (GeneFromGTF g) {
		List<GTFRecord> result = new ArrayList<GTFRecord> ();
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
                        CONSENSUS_INTRON_FEATURE_TYPE));
            }
        }
        records.addAll(consensusIntrons);
    }

	public GTFRecord getGTFRecord(GeneFromGTF g) {
		GTFRecord r = new GTFRecord(g.getSequence(), g.getStart(), g.getEnd(), g.isNegativeStrand(), 
				g.getGeneID(), g.getName(), null, null, g.getTranscriptType(), GENE_FEATURE_TYPE);
		return (r);
	}
	
	
	public List<GTFRecord> enhanceTranscript (GeneFromGTF.TranscriptFromGTF t) {
		List<GTFRecord> result = new ArrayList<GTFRecord> ();
		GeneFromGTF g = t.getGene();
		GTFRecord transcriptRecord = new GTFRecord(g.getSequence(), t.codingStart, t.codingEnd, g.isNegativeStrand(), g.getGeneID(), g.getName(), t.getTranscriptName(), t.getTranscriptID(), t.getTranscriptType(), TRANSCRIPT_FEATURE_TYPE);
		result.add(transcriptRecord);
		
		List<Interval> introns = getIntronIntervals(getIntervals(t.exons));
		List<GTFRecord> intronR = getGTFRecordsFromIntronIntervals(introns, t);
		result.addAll(intronR);
		List<GTFRecord> exons = getGTFRecordsFromExons(t);
		result.addAll(exons);
		return (result);
		
	}
	
	
	List<GTFRecord> getGTFRecordsFromExons(GeneFromGTF.TranscriptFromGTF t) {
		List<GTFRecord> result = new ArrayList<GTFRecord> ();
		GeneFromGTF g = t.getGene();
		for (Exon e: t.exons) {
			GTFRecord exonRecord = new GTFRecord(g.getSequence(), e.start, e.end, g.isNegativeStrand(), g.getGeneID(), g.getName(), t.getTranscriptName(), t.getTranscriptID(), t.getTranscriptType(), EXON_FEATURE_TYPE);
			result.add(exonRecord);
		}
		return (result);
	}
	
	List<GTFRecord> getGTFRecordsFromIntronIntervals (List<Interval> introns, GeneFromGTF.TranscriptFromGTF t) {
		List<GTFRecord> result = new ArrayList<GTFRecord> ();
		GeneFromGTF g = t.getGene();
		for (Interval i: introns) {
			GTFRecord intronRecord = new GTFRecord(g.getSequence(), i.getStart(), i.getEnd(), g.isNegativeStrand(), g.getGeneID(), g.getName(), t.getTranscriptName(), t.getTranscriptID(), t.getTranscriptType(), INTRON_FEATURE_TYPE);
			result.add(intronRecord);
		}		
		return (result);
	}
	
	
	List<Interval> getIntronIntervals(List<Interval> exons) {
		if (exons.size()==1) return new ArrayList<Interval>(0);
		
		List<Interval> result = new ArrayList<Interval>(exons.size()-1);
		for (int i=0; i<exons.size()-1; i++) {
			int s = exons.get(i).getEnd()+1;
			int e = exons.get(i+1).getStart()-1;
			Interval intron = new Interval(null, s, e);
			result.add(intron);
		}
		return (result);
	}
	
	List<Interval> getIntervals (Exon [] exons) {
		List<Interval> result = new ArrayList<Interval>();
		for (Exon e: exons) {
			result.add(new Interval(null, e.start, e.end));
		}
		return result;
	}
	
	OverlapDetector<GeneFromGTF> getOverlapDetector(List<GTFRecord> records) {
		Map<String, List<GTFRecord>> recordsByGene = getRecordsByGene(records); 
		// this isn't pretty, think about a better pattern for exposing this functionality in GTFReader...
		List<SAMSequenceRecord> fake=  new ArrayList<SAMSequenceRecord>();
		GTFReader r = new GTFReader(new File(""), new SAMSequenceDictionary(fake));
		OverlapDetector<GeneFromGTF> od = r.convert(recordsByGene);
		return (od);
	}
	
	
	Map<String, List<GTFRecord>> getRecordsByGene (List<GTFRecord> records) {
		Map<String, List<GTFRecord>> recordsByGene=new HashMap<String, List<GTFRecord>>();
		for (GTFRecord r: records) {
			List<GTFRecord> l = recordsByGene.get(r.getGeneID());
			if (l==null) {
				l = new ArrayList<GTFRecord>();
				recordsByGene.put(r.getGeneID(),l);
			}
			l.add(r);
		}
		return (recordsByGene);
		
	}
	
}
