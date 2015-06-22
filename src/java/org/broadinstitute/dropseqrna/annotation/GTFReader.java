package org.broadinstitute.dropseqrna.annotation;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;

import java.io.File;
import java.util.*;

import org.broadinstitute.dropseqrna.annotation.GeneFromGTF.TranscriptFromGTF;

import picard.annotation.AnnotationException;
import picard.util.TabbedTextFileWithHeaderParser;

/**
 * Loads gene annotations from a GTF file into an OverlapDetector<Gene>.  Discards annotations that are not
 * internally consistent, e.g. transcripts on different chromosomes or different strands.
 * This borrows heavily from RefFlatReader.  Thanks picard!
 * The big difference in GTF vs RefFlat is that refFlat defines a transcript on a single line with all exons, while GTF gives you a 1 exon per line with the gene/transcript it's associated with.
 * This forces you to read 1 or more lines to assemble a transcript properly.  Coding start/stop information is also optional, and if it exists also exists on separate lines.
 * @see http://www.sanger.ac.uk/resources/software/gff/spec.html#t_2
 */
public class GTFReader {
	
    private static final Log LOG = Log.getInstance(GTFReader.class);
    
    
    // These are in the order that columns appear in refFlat format.
    public enum GTFFlatColumns{CHROMOSOME, SOURCE, FEATURE, START, END, SCORE, STRAND, FRAME, ATTRIBUTE}
    
    private static final String[] GTFColumnLabels = new String[GTFFlatColumns.values().length];
    
    static {
        for (int i = 0; i < GTFColumnLabels.length; ++i) {
            GTFColumnLabels[i] = GTFFlatColumns.values()[i].name();
        }
    }

    private final File gtfFlatFile;
    private final SAMSequenceDictionary sequenceDictionary;
    // Keep track of errors already reported to reduce verbosity
    private final Set<String> skippedChromosomeTranscriptDescription = new HashSet<>();
    private final Set<String> unrecognizedSequences = new HashSet<>();

    public GTFReader(final File gtfFlatFile, final SAMSequenceDictionary sequenceDictionary) {
        this.gtfFlatFile = gtfFlatFile;
        this.sequenceDictionary = sequenceDictionary;
    }
    
    public GTFReader (final File gtfFlatFile, final File sequenceDictionary) {
    	final SAMFileReader r = new SAMFileReader(sequenceDictionary);
    	SAMSequenceDictionary dict = r.getFileHeader().getSequenceDictionary();
    	r.close();
    	this.gtfFlatFile=gtfFlatFile;
    	this.sequenceDictionary=dict;
    }
    
    
    
    static OverlapDetector<GeneFromGTF> load(final File refFlatFile, final SAMSequenceDictionary sequenceDictionary) {
        return new GTFReader(refFlatFile, sequenceDictionary).load();
    }
	
    public OverlapDetector<GeneFromGTF> load() {
        final TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(gtfFlatFile, GTFColumnLabels);
        Map<String, List<GTFRecord>> gtfLinesByGene = getRecordsByGene(parser);
        final OverlapDetector<GeneFromGTF> overlapDetector =  convert(gtfLinesByGene);
        return overlapDetector;
    }
    
    public OverlapDetector<GeneFromGTF> convert(Map<String, List<GTFRecord>> gtfLinesByGene) {
    	final OverlapDetector<GeneFromGTF> overlapDetector = new OverlapDetector<GeneFromGTF>(0, 0);
    	int longestInterval = 0;
        int numIntervalsOver1MB = 0;

        for (final List<GTFRecord> transcriptLines : gtfLinesByGene.values()) {
            try {
                final GeneFromGTF gene = makeGeneFromGTFLines(transcriptLines);
                overlapDetector.addLhs(gene, gene);
                if (gene.length() > longestInterval) longestInterval = gene.length();
                if (gene.length() > 1000000) ++numIntervalsOver1MB;
            } catch (AnnotationException e) {
                LOG.info(e.getMessage() + " -- skipping");
            }
        }
        LOG.debug("Longest gene: " + longestInterval + "; number of genes > 1MB: " + numIntervalsOver1MB);
        LOG.debug("Total number of genes loaded [" + overlapDetector.getAll().size() +"]");
        return overlapDetector;
    }

    public Set<String> getUnrecognizedSequences() {
        return unrecognizedSequences;
    }

    private boolean isSequenceRecognized(final String sequence) {
        return (sequenceDictionary.getSequence(sequence) != null);
    }


    private GeneFromGTF makeGeneFromGTFLines(final List<GTFRecord> transcriptLines) {
    	GTFRecord lineOne=transcriptLines.get(0);
    	
        String geneName=lineOne.getGeneName();
       
        final boolean transcriptNegStrand= lineOne.isNegativeStrand();

        // Figure out the extend of the gene
        int start = Integer.MAX_VALUE;
        int end = Integer.MIN_VALUE;
        for (final GTFRecord r: transcriptLines) {
            start = Math.min(start, r.getStart());
            end   = Math.max(end,   r.getEnd());
        }

        final GeneFromGTF gene = new GeneFromGTF(lineOne.getChromosome(), start, end, transcriptNegStrand, geneName, lineOne.getFeatureType(), 
        		lineOne.getGeneID(), lineOne.getTranscriptType());
        final Map<String, List<GTFRecord>> gtfLinesByTranscript = getLinesByTranscript(transcriptLines, gene);
        for (String k: gtfLinesByTranscript.keySet()) {
        	List<GTFRecord> tl = gtfLinesByTranscript.get(k);
        	addTranscriptToGeneFromGTFLines(gene, tl);
        }
        
        
        return gene;
    }

    /**
     * Conversion from 0-based half-open to 1-based inclusive intervals is done here.
     */
    private TranscriptFromGTF addTranscriptToGeneFromGTFLines(final GeneFromGTF gene, final List<GTFRecord> transcriptLines) {
    	
        final String geneName = gene.getName();
        
        GTFRecord lineOne=transcriptLines.get(0);
    	final String transcriptName = lineOne.getTranscriptName();
    	final String transcriptID = lineOne.getTranscriptID();
        final String transcriptDescription = geneName + ":" + transcriptName;
        String transcriptType = lineOne.getTranscriptType();
        
        List<Exon> exons = new ArrayList<Exon>();
        
        int transcriptionStart = Integer.MAX_VALUE;
        int transcriptionEnd = Integer.MIN_VALUE;
        int codingStart = Integer.MAX_VALUE;
        int codingEnd = Integer.MIN_VALUE;
        
        for (GTFRecord r: transcriptLines) {
        	String featureType=r.getFeatureType();
        	int start = r.getStart();
    		int end = r.getEnd();
    		
        	if (featureType.equals("exon")) {
        		Exon e = new Exon(start,end);
                exons.add(e);
                transcriptionStart = Math.min(transcriptionStart, start);
                transcriptionEnd = Math.max(transcriptionEnd, end);

            }
            if (featureType.equals("CDS")) {
            	codingStart=Math.min(codingStart, start);
        		codingEnd=Math.max(codingEnd, end);
            }           
        }
        
        if (codingStart==Integer.MAX_VALUE) codingStart=transcriptionStart;
        if (codingEnd==Integer.MIN_VALUE) codingEnd=transcriptionEnd;
        
        final TranscriptFromGTF tx = gene.addTranscript(transcriptName, transcriptionStart, transcriptionEnd, codingStart, codingEnd, exons.size(), transcriptName, transcriptID, transcriptType);
        
        for (int i = 0; i < exons.size(); ++i) {
        	Exon e = exons.get(i);
            if (e.start > e.end) {
                throw new AnnotationException("Exon has 0 or negative extent for " + transcriptDescription);
            }
            if (i > 0 && exons.get(i-1).end >= exons.get(i).start) {
                throw new AnnotationException("Exons overlap for " + transcriptDescription);
            }
            tx.addExon(e.start, e.end);
        }
        
   
        return tx;
        
    	
    }
    
    private Map<String, List<GTFRecord>> getRecordsByGene (TabbedTextFileWithHeaderParser parser) {
    	final int expectedColumns = GTFFlatColumns.values().length;
    	
    	final Map<String, List<GTFRecord>> gtfRecordsByGene =new HashMap<String, List<GTFRecord>>();

        for (final TabbedTextFileWithHeaderParser.Row row : parser) {
            final int lineNumber = parser.getCurrentLineNumber(); // getCurrentLineNumber returns the number of the next line
            if (lineNumber%100000==0) {
            	LOG.info("Processed " + lineNumber +  " lines");
            }
            if (row.getFields().length != expectedColumns) {
                throw new AnnotationException("Wrong number of fields in GTF file " + gtfFlatFile + " at line " +
                        lineNumber);
            }
            
            GTFRecord r = parseLine(row);
            final String transcriptDescription = r.getGeneName() + ":" + r.getTranscriptName();
            
            if (!isSequenceRecognized(r.getChromosome())) {
                unrecognizedSequences.add(r.getChromosome());
                if (skippedChromosomeTranscriptDescription.add(transcriptDescription + "\t" + r.getChromosome())) {
                    LOG.debug("Skipping " + transcriptDescription + " due to unrecognized sequence " + r.getChromosome());
                }
            } else {
                List<GTFRecord> transcriptRecords = gtfRecordsByGene.get(r.getGeneName());
                if (transcriptRecords == null) {
                	transcriptRecords = new ArrayList<GTFRecord>();
                    gtfRecordsByGene.put(r.getGeneName(), transcriptRecords);
                }
                transcriptRecords.add(r);
            }
        }
        return gtfRecordsByGene;
    }
    
    private GTFRecord parseLine (TabbedTextFileWithHeaderParser.Row row) {
    	GTFRecord r = null;
    	String attributes= row.getField(GTFFlatColumns.ATTRIBUTE.name());
        Map<String, String> attributesMap = AnnotationUtils.getInstance().parseOptionalFields(attributes);
        
        String chromosome = row.getField(GTFFlatColumns.CHROMOSOME.name());
        int start = Integer.parseInt(row.getField(GTFFlatColumns.START.name()));
        int end = Integer.parseInt(row.getField(GTFFlatColumns.END.name()));
        String strand = row.getField(GTFFlatColumns.STRAND.name());
        String featureType = row.getField(GTFFlatColumns.FEATURE.name());
        
        String geneName=attributesMap.get("gene_name");
        String geneID=attributesMap.get("gene_id");
		String transcriptName=attributesMap.get("transcript_name");
		String transcriptID=attributesMap.get("transcript_id");
		String transcriptType=attributesMap.get("gene_biotype"); 
			
		boolean negativeStrand=false;
		if (strand.equals("-")) negativeStrand=true;
		
        r=new GTFRecord (chromosome, start, end, negativeStrand, geneID, geneName, transcriptName, transcriptID, transcriptType, featureType);
    	return (r);
    }
    
    
    
    private Map<String, List<GTFRecord>> getLinesByTranscript (List<GTFRecord> transcriptLines, GeneFromGTF gene) {
    	String chromosome = gene.getSequence();
    	boolean geneIsPositiveStrand = gene.isPositiveStrand();
    	String geneName = gene.getName();
    	
    	final Map<String, List<GTFRecord>> gtfLinesByTranscript = new HashMap<String, List<GTFRecord>>();
        
        for (final GTFRecord record: transcriptLines) {
            if (geneIsPositiveStrand==record.isNegativeStrand()) {
                throw new AnnotationException("Strand disagreement in GTF file for gene " + geneName);
            }
            if (!gene.getSequence().equals(record.getChromosome())) {
                throw new AnnotationException("Chromosome disagreement(" + chromosome + " != " + record.getChromosome() +
                                                      ") in GTF file for gene " + geneName);
            }
            if (record.getFeatureType().equals("gene")) {
                if (record.getChromosome().equals(gene.getContig()) &&
                        record.getStart() == gene.getStart() &&
                        record.getEnd() == gene.getEnd() &&
                        record.isNegativeStrand() == gene.isNegativeStrand()) {
                    // Do not include this in lines by transcript.
                    continue;
                } else {
                    throw new RuntimeException(String.format("gene GTFRecord(%s) != GeneFromGTF(%s)", record.toString(), gene.toString()));
                }
            }
            
            String transcriptID = record.getTranscriptID();
            if (transcriptID == null) {
                throw new RuntimeException("GTFRecord does not have transcriptID: " + record);
            }
            //parse out the transcripts by ID and group them.
            List<GTFRecord> tl = gtfLinesByTranscript.get(transcriptID);
            if (tl == null) {
                tl = new ArrayList<GTFRecord>();
                gtfLinesByTranscript.put(transcriptID, tl);
            }
            tl.add(record);
        }
        
        // sort by genomic position.
        for (String k: gtfLinesByTranscript.keySet()) {
        	List<GTFRecord> v = gtfLinesByTranscript.get(k);
        	Collections.sort(v);
        	v.size();
        }
        
    	return (gtfLinesByTranscript);
    }
    
    public class Exon {
        public final int start;
        public final int end;

        public Exon(final int start, final int end) {
            this.start = start;
            this.end = end;
        }
    }
}
