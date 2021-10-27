package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import picard.util.TabbedTextFileWithHeaderParser;
import picard.util.TabbedTextFileWithHeaderParser.Row;

/**
 * Reads GatherDigitalAlleleCounts allele frequency output into two different forms:
 * 1) A simple map from the interval to the UMI minor allele frequency
 * 2) Reconstruct the list of GdacAlleleFrequency objects that generated this file 
 * @author nemesh
 *
 */
public class GdacAlleleFrequencyReader {

	private final File inputFile;
	private static final Log log = Log.getInstance(GdacAlleleFrequencyReader.class);
	private List<String> expectedHeaderFull = new ArrayList<String>(Arrays.asList("chromosome", "position", "ref_allele", "alt_allele", 
			"ref_reads", "alt_reads", "ref_umi", "alt_umi", "maf_reads", "maf_umi"));
	
	private List<String> expectedHeaderMinimum = new ArrayList<String>(Arrays.asList("chromosome", "position", "maf_umi"));
	
	
	public GdacAlleleFrequencyReader (final File inputFile) {
		this.inputFile=inputFile;	
	}
	
	/**
	 * Simple map of Interval to minor allele frequency (as determined by UMI counts)
	 * @return Map of Interval to minor allele frequency (as determined by UMI counts)
	 */
	public Map<Interval, Double> getUmiAlleleFrequencyMap () {
		TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(this.inputFile);
		validateHeader(parser.getColumnNames(), expectedHeaderMinimum);
		
		Map<Interval, Double> result = new HashMap<>();
		
		Iterator<Row> iter = parser.iterator();
		
 		// parse each line into the map
		while (iter.hasNext()) {
			Row r = iter.next();
			Interval i = parseInterval(r);
			double maf = Double.parseDouble(r.getField("maf_umi"));
			result.put(i, maf);						
		}
		parser.close();		
		return result;
	}
	
	private Interval parseInterval(Row r) {
		int pos = Integer.parseInt(r.getField("position"));
		Interval i = new Interval (r.getField("chromosome"), pos, pos);
		return (i);
	}
	
	/**
	 * Parse the input file, and return a list of GdacAlleleFrequency objects
	 * @return The GdacAlleleFrequency list that generated the input file.
	 */
	public List<GdacAlleleFrequency> parseFile () {
		TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(this.inputFile);
		validateHeader(parser.getColumnNames(), expectedHeaderFull);
		
		List<GdacAlleleFrequency> result = new ArrayList<>();
		
		Iterator<Row> iter = parser.iterator();
		
 		// parse each line into the map
		while (iter.hasNext()) {
			Row r = iter.next();
			Interval snpInterval = parseInterval(r);
						
			char refAllele=r.getField("ref_allele").charAt(0);
			char altAllele=r.getField("alt_allele").charAt(0);
			
			int refReadCount = Integer.parseInt(r.getField("ref_reads"));
			int altReadCount = Integer.parseInt(r.getField("alt_reads"));
			int refUmiCount = Integer.parseInt(r.getField("ref_umi"));
			int altUmiCount = Integer.parseInt(r.getField("alt_umi"));
			GdacAlleleFrequency freq=new GdacAlleleFrequency(snpInterval, refAllele, altAllele, refReadCount, altReadCount, refUmiCount, altUmiCount);
			result.add(freq);						
		}		
		parser.close();
		return result;
	}
	
	private void validateHeader (Set<String> header, List<String> expectedHeader) {
		// validate that the expected header lines are contained				
		if (!header.containsAll(expectedHeader)) {			
			log.error("Parsing cell contamination map, expected header values" + expectedHeader.toString());
			log.error("File header ", header.toString());
			throw new IllegalArgumentException("Cell contamination file has the wrong headers!");
		}
	}
}
