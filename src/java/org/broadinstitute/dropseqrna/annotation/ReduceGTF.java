package org.broadinstitute.dropseqrna.annotation;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.dropseqrna.cmdline.MetaData;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.TabbedInputParser;

/**
 * GTF files are annoyingly complex with a poor definition of what data is in them.
 * So hey, why not write a file parser.
 * This program reduces the GTF file in to a simplier, easier to parse format, while simultaneously allowing for data to be filtered.
 * @author nemesh
 *
 */

@CommandLineProgramProperties(
        usage = "GTF files are annoyingly complex with a poor definition of what data is in them. So hey, why not write a file parser. This program reduces the GTF file in to a simplier, easier to parse format, while simultaneously allowing for data to be filtered.",
        usageShort = "Parse and simplify a GTF file into an easier to use format.",
        programGroup = MetaData.class
)
public class ReduceGTF extends CommandLineProgram {

    private static final Log log = Log.getInstance(ReduceGTF.class);
    private static final List<String> DEFAULT_FEATURE_TYPES = CollectionUtil.makeList("gene", "transcript", "exon");
    private static final List<String> DEFAULT_IGNORED_FUNC_TYPES = CollectionUtil.makeList(
            "pseudogene", "polymorphic_pseudogene", "TR_J_pseudogene", "TR_V_pseudogene", "IG_C_pseudogene",
            "IG_J_pseudogene", "IG_V_pseudogene");
    private static final String NA = "NA";

    @Option(shortName = StandardOptionDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, doc="The reference sequence dictionary." +
            "  Only chromosomes found in this file AND the GTF file will be retained.")
	public File SEQUENCE_DICTIONARY;
	
	@Option(doc="The GTF file to reduce")
	public File GTF;
	
	@Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME,doc="The output reduced GTF file.")
	public File OUTPUT;
	
	@Option(doc="The field name that contains the function of the gene/transcript")
	public String FUNCTION_FIELD="gene_biotype";
	
	@Option(doc="Feature type(s) to extract. Only lines of the GTF that have these feature types will be extracted.  " +
            "This is the 3rd field of the GTF file, some examples of standard feature types are CDS, start_codon, stop_codon, and exon. ")
	public List<String> FEATURE_TYPE = DEFAULT_FEATURE_TYPES;
	
	@Option(doc="Functional type(s) to ignore.  These are values in the FUNCTIONAL_FIELD column in the GTF file.")
	public List<String> IGNORE_FUNC_TYPE = DEFAULT_IGNORED_FUNC_TYPES;
	
	@Option(doc="Enhance this reduced GTF file with genes,transcripts,introns, and consensus introns.  This is real " +
            "handy when your GTF file only defines exons, but has the transcript and gene IDs they belong to.")
	public boolean ENHANCE_GTF=true;
	
	@Override
	protected int doWork() {
		IOUtil.assertFileIsReadable(SEQUENCE_DICTIONARY);
		IOUtil.assertFileIsReadable(GTF);
		IOUtil.assertFileIsWritable(this.OUTPUT);
		
		Set<String> featureTypes = new HashSet<>(FEATURE_TYPE);
		Set<String> ignoredFunctionalTypes = new HashSet<>(IGNORE_FUNC_TYPE);
		
		List<GTFRecord> records = parseGTF (this.GTF, this.SEQUENCE_DICTIONARY, featureTypes, ignoredFunctionalTypes, true );
		
		PrintStream out = new PrintStream(IOUtil.openFileForWriting(OUTPUT));
		writeHeader(out);
		
		// if no enhancement needed, just write out the results.
		if (!ENHANCE_GTF) {	
			writeRecords (out, records);
			out.close();
			return 0;
		}
		
		// enhance
		EnhanceGTFRecords e = new EnhanceGTFRecords();
		records = e.enhanceGTFRecords(records);
		
		//sort the data before writing?		
		Collections.sort(records, new GenomicOrderComparator(getDict(this.SEQUENCE_DICTIONARY)));
		
		// write results
		writeRecords(out, records);
		
		out.close();
		return 0;
	}

	
	
	
	
	
	List<GTFRecord> parseGTF (File gtfFile, File dictFile, Set<String> featureTypes, Set<String> ignoredFunctionalTypes, boolean reportProgress ) {
		SAMSequenceDictionary dict = getDict(dictFile);
		TabbedInputParser parser = new TabbedInputParser(false, gtfFile);
		List<GTFRecord> records=new ArrayList<GTFRecord>();
		
		int counter=0;
		while (parser.hasNext()) {
			if (counter%100000==0 && reportProgress) log.info("Progress [" + counter +"] lines");
			String [] line = parser.next();
			
			GTFRecord r = parseLine(line, this.FUNCTION_FIELD, featureTypes, ignoredFunctionalTypes, dict);
			
			if (r!=null) {	
				records.add(r);	
			}
							
			counter++;
		}
		parser.close();
		return (records);
	}
	
	
	private GTFRecord parseLine (String [] line, String functionField, Set<String> featureTypes, Set<String> ignoredFunctionalTypes, SAMSequenceDictionary dict) {
		String chr = line[0];
		SAMSequenceRecord index = dict.getSequence(chr);
		if (index==null) {
			return (null);
		}
		
		String featureType=line[2];
		// throw away lines that don't have the right feature type
		if (featureTypes.contains(featureType)==false) {
			return (null);
		}
		
		int startPos=Integer.parseInt(line[3]);
		int endPos=Integer.parseInt(line[4]);
		String strand=line[6];
		
		Map<String, String> optionalAnnosMap= parseOptionalFields(line[8]);
		
		String geneName=optionalAnnosMap.get("gene_name");
		String transcriptName=optionalAnnosMap.get("transcript_name");
		String transcriptType=optionalAnnosMap.get(functionField);
		String geneID=optionalAnnosMap.get("gene_id");
		String transcriptID=optionalAnnosMap.get("transcript_id");
		if (ignoredFunctionalTypes.contains(transcriptType)) {
			return (null);
		}
		
		boolean negativeStrand=false;
		if (strand.equals("-")) negativeStrand=true;
		
		GTFRecord r = new GTFRecord(chr, startPos, endPos, negativeStrand, geneID, geneName, transcriptName, transcriptID, transcriptType, featureType);
		
		return (r);
	}
	
	private void writeHeader (PrintStream out) {
		String [] line = {"chr", "start", "end", "strand", "gene_name", "gene_id", "transcript_name", "transcript_id",
                "transcriptType", "annotationType"};
		String h = StringUtils.join(line, "\t");
		out.println(h);
	}
	
	private void writeRecords (PrintStream out, List<GTFRecord> records) {
		for (GTFRecord r: records) {
			writeLine(r, out);
		}
	}
	
	Map<String, String> parseOptionalFields(String optional) {
		Map<String, String> result = new HashMap<String, String>();
		String [] o = optional.split(";");
		for (String s: o) {
			s=s.replaceAll("\"", "");
			s=s.trim();
			String [] z= s.split(" ");
			String k = z[0];
			String v = z[1];
			result.put(k, v);
		}
		return (result);
	}
	
	private void writeLine (GTFRecord r, PrintStream out) {
		if (r==null) return;
		String [] line={r.getChromosome(),new Integer(r.getStart()).toString(), new Integer(r.getEnd()).toString(), r.getStrandAsString(), r.getGeneName(), r.getGeneID(),
				r.getTranscriptName(), r.getTranscriptID(), r.getTranscriptType(), r.getFeatureType()};
        for (int i = 0; i < line.length; ++i) {
            if (line[i] == null || line[i].isEmpty()) {
                line[i] = NA;
            }
        }
		String h = StringUtils.join(line, "\t");
		out.println(h);
	}
	
	private SAMSequenceDictionary getDict (File file) {
		final SAMFileReader r = new SAMFileReader(file);
		SAMSequenceDictionary dict = r.getFileHeader().getSequenceDictionary();
		r.close();
		return (dict);
	}
	
	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new ReduceGTF().instanceMain(args));
	}
	
	
}
