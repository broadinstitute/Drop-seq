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

import java.io.File;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.OverlapDetector;
import picard.annotation.Gene;

/**
 * This class handles reading both GTFs and refFlat files.
 * @author nemesh
 *
 */
public class GeneAnnotationReader {

	public static OverlapDetector<Gene> loadAnnotationsFile(final File annotationFile, final SAMSequenceDictionary sequenceDictionary) {

		// trim off the potential compression tags on the file.
		String f = annotationFile.getName();
		if (f.endsWith(".bz2")) f=f.replaceAll("\\.bz2", "");
		if (f.endsWith(".gz")) f=f.replaceAll("\\.gz", "");

		if (f.endsWith(".gtf"))
			return loadGTFFile(annotationFile, sequenceDictionary);
		else if (f.endsWith(".refFlat"))
			return loadRefFlat(annotationFile, sequenceDictionary);

		throw new SAMException("Unable to determine file format for gene annotation file.  Expect [.gtf | .refFlat]");
	}


	public static OverlapDetector<Gene> loadAnnotationsFile(final File annotationFile, final File sequenceDictionary) {
		SAMSequenceDictionary sd = getSequenceDictionary(sequenceDictionary);
		return (loadAnnotationsFile(annotationFile, sd));
	}


	private static OverlapDetector<Gene> loadGTFFile(final File gtfFile, final SAMSequenceDictionary sequenceDictionary) {
		@SuppressWarnings({ "rawtypes", "unchecked" })
		OverlapDetector<Gene> result = (OverlapDetector) GTFReader.load(gtfFile, sequenceDictionary);
		return result;
	}

	private static OverlapDetector<Gene> loadRefFlat(final File refFlatFile, final SAMSequenceDictionary sequenceDictionary) {
		return picard.annotation.GeneAnnotationReader.loadRefFlat(refFlatFile, sequenceDictionary);
    }

	private static SAMSequenceDictionary getSequenceDictionary (final File sd) {
		SamReader reader = SamReaderFactory.makeDefault().open(sd);
		SAMSequenceDictionary dict = reader.getFileHeader().getSequenceDictionary();
		CloserUtil.close(reader);
    	return dict;
	}
}
