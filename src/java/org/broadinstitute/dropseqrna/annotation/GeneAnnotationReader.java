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
