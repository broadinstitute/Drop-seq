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

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.dropseqrna.cmdline.MetaData;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import picard.annotation.Gene;
import picard.annotation.Gene.Transcript.Exon;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.PrintStream;


/**
 * Converts an annotations file to refFlat format.
 * Yes, you can read in RefFlat and output the same file.  I dunno why, but go nuts.
 * @author nemesh
 *
 */

@CommandLineProgramProperties(
        usage = "If you really want to have refFlat files instead of GTF files, then this is for you.  Pretty handy for Picard tools that accept refFlat instead of GTF files.",
        usageShort = "Convert various annotations formats to RefFlat",
        programGroup = MetaData.class
)

public class ConvertToRefFlat extends CommandLineProgram {

	
	@Option(doc="The annotations set to use to label the read.  This can be a GTF or a refFlat file.")
	public File ANNOTATIONS_FILE;
	
	@Option(shortName = StandardOptionDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, doc="The reference sequence dictionary.  Only chromosomes found in this file AND the annotations file will be retained.")
	public File SEQUENCE_DICTIONARY;
	
	@Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output refFlat file")
	public File OUTPUT;
	
	@Override
	protected int doWork() {
		IOUtil.assertFileIsReadable(this.ANNOTATIONS_FILE);
		IOUtil.assertFileIsWritable(this.OUTPUT);
		
		OverlapDetector<Gene> od = GeneAnnotationReader.loadAnnotationsFile(this.ANNOTATIONS_FILE, this.SEQUENCE_DICTIONARY);
		
		PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(OUTPUT));
		
		for (Gene g: od.getAll()) {
			for (Gene.Transcript t: g) {
				RefFlatRecord r = convertTranscript(t);
				out.println(r.toString());
			}
		}
		out.close();
		return 0;
	}
	
	
	
	RefFlatRecord convertTranscript (Gene.Transcript t) {
		Gene g = t.getGene();
		
		RefFlatRecord r = new RefFlatRecord(g.getName(), t.name, g.getContig(), g.isNegativeStrand(), t.transcriptionStart,
				t.transcriptionEnd, t.codingStart, t.codingEnd);
		
		for (Exon e: t.exons) {
			r.addExonStart(e.start);
			r.addExonEnd(e.end);
		}
		return (r);
	}
	
	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new ConvertToRefFlat().instanceMain(args));
	}
	
	
}
