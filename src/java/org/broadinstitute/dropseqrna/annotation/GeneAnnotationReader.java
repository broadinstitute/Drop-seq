package org.broadinstitute.dropseqrna.annotation;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.OverlapDetector;

import java.io.File;

import picard.annotation.Gene;

/**
 * This class handles reading both GTFs and refFlat files.
 * @author nemesh
 *
 */
public class GeneAnnotationReader {

	public static OverlapDetector<Gene> loadAnnotationsFile(File annotationFile, SAMSequenceDictionary sequenceDictionary) {
		
		// trim off the potential compression tags on the file.
		String f = annotationFile.getName();
		if (f.endsWith(".bz2")) f=f.replaceAll("\\.bz2", "");
		if (f.endsWith(".gz")) f=f.replaceAll("\\.gz", "");
		
		if (f.endsWith(".gtf")) {
			return loadGTFFile(annotationFile, sequenceDictionary);
		} else if (f.endsWith(".refFlat")) {
			return loadRefFlat(annotationFile, sequenceDictionary);
		} 
		
		throw new SAMException("Unable to determine file format for gene annotation file.  Expect [.gtf | .refFlat]");        
	}

	
	public static OverlapDetector<Gene> loadAnnotationsFile(File annotationFile, File sequenceDictionary) {
		SAMSequenceDictionary sd = getSequenceDictionary(sequenceDictionary);
		return (loadAnnotationsFile(annotationFile, sd));
	}
	
	
	private static OverlapDetector<Gene> loadGTFFile(File gtfFile, SAMSequenceDictionary sequenceDictionary) {
		@SuppressWarnings({ "rawtypes", "unchecked" })
		OverlapDetector<Gene> result = (OverlapDetector) GTFReader.load(gtfFile, sequenceDictionary);
		return result;
	}
	
	private static OverlapDetector<Gene> loadRefFlat(File refFlatFile, SAMSequenceDictionary sequenceDictionary) {
		return picard.annotation.GeneAnnotationReader.loadRefFlat(refFlatFile, sequenceDictionary);
    }
	
	private static SAMSequenceDictionary getSequenceDictionary (File sd) {
		final SAMFileReader r = new SAMFileReader(sd);
    	SAMSequenceDictionary dict = r.getFileHeader().getSequenceDictionary();
    	r.close();
    	return dict;
	}
}
