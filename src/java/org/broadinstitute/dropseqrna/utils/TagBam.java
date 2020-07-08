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
package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.*;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.util.Iterator;

@CommandLineProgramProperties(summary = "Adds a BAM tag to every read",
        oneLineSummary = "Adds a BAM tag to every read",
        programGroup = DropSeq.class)
public class TagBam extends CommandLineProgram{
	private final Log log = Log.getInstance(TagBam.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.")
	public File INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output bam")
	public File OUTPUT;
	
	@Argument(doc="Tag name to set on every read")
	public String TAG_NAME;
	
	@Argument(doc="Tag value to set on every read.  If this value is not set, the TAG is removed from all records of the BAM.", optional=true)
	public String TAG_VALUE;
	
	@Override
	protected int doWork() {
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		SamReader reader = SamReaderFactory.makeDefault().open(INPUT);
		
		SAMFileHeader h= reader.getFileHeader();
		SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(h, true, OUTPUT);
		Iterator<SAMRecord> iter = reader.iterator();
		ProgressLogger p = new ProgressLogger(this.log);
		
		while (iter.hasNext()) {
			SAMRecord r = iter.next();
			r.setAttribute(TAG_NAME, TAG_VALUE);
			writer.addAlignment(r);
			p.record(r);
		}
		
		CloserUtil.close(reader);
		writer.close();
		
		return (0);
	}
	
	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new TagBam().instanceMain(args));
	}
}
