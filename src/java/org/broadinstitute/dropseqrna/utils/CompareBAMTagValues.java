/*
 * MIT License
 *
 * Copyright 2018 Broad Institute
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

import java.io.File;
import java.util.Iterator;
import java.util.List;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Tests that two BAMs have the same TAG values per read.  "
		+ "We assume the two BAMs tested have the same set of reads in the same order, and test the read names "
		+ "for each iteration match, and the values for the TAGs listed are equal.",
oneLineSummary = "Tests that two BAMs have the same TAG values per read.",
programGroup = DropSeq.class)

public class CompareBAMTagValues extends CommandLineProgram {

	private final Log log = Log.getInstance(CompareBAMTagValues.class);

	@Argument(doc = "The input SAM or BAM file to analyze.", optional=false)
	public File INPUT_1;

	@Argument(doc = "The input SAM or BAM file to analyze.", optional=false)
	public File INPUT_2;

	@Argument(doc="A list of tags to test.  Values of these tags will be compared across BAMs for equality.", optional=false)
	public List<String> TAGS;

	@Override
	public int doWork() {
		IOUtil.assertFileIsReadable(INPUT_1);
		IOUtil.assertFileIsReadable(INPUT_2);

		Iterator<SAMRecord> iterator1 = SamReaderFactory.makeDefault().open(INPUT_1).iterator();
		Iterator<SAMRecord> iterator2 = SamReaderFactory.makeDefault().open(INPUT_2).iterator();

		while (iterator1.hasNext() && iterator2.hasNext()) {
			SAMRecord r1 = iterator1.next();
			SAMRecord r2 = iterator2.next();
			// log.info(r1.toString() + " tags read 1: " + getTagsAsString(r1, this.TAGS));
			boolean test = compareTags(r1, r2, this.TAGS) ;
			if (!test) {
				log.info("Difference found for read " + r1.toString()+ " read 1: " + getTagsAsString(r1, this.TAGS)+ " read 2: "+ getTagsAsString(r2, this.TAGS));				
				return 1;
			}
		}
		log.info("DONE");
		return 0;
	}
	
	private String getTagsAsString (SAMRecord r, List<String> tags) {
		StringBuilder b = new StringBuilder();
		for (String t: tags) {
			b.append(t + "[" +  r.getStringAttribute(t) +"] ");
		}
		return b.toString();
	}

	boolean compareTags (final SAMRecord r1, final SAMRecord r2, final List<String> tags) {
		boolean sameRead = r1.getReadName().equals(r2.getReadName());
		if (!sameRead)
			log.error("Read names differ at iteration.  R1: "+ r1.toString(), "R2: ", r2.toString());

		for (String tag: tags) {
			Object o1 = r1.getAttribute(tag);
			Object o2 = r2.getAttribute(tag);
			if (o1==null && o2==null) return true;
			if ((o1==null && o2!=null ) || (o1!=null && o2==null)) {
				log.error("Read tag values differ for tag: ["+ tag.toString()+ "] R1 is null [" + String.valueOf(o1==null) + "] R2 is null [" + String.valueOf(o2==null)+"]");
				return false;
			}


			if (!o1.equals(o2)) {
				log.error("Read tag values differ for tag: ["+ tag.toString()+ "] "  + r1.toString()+ "[" +o1.toString()+ "] "+ r2.toString() +" [" + o2.toString()+"]");
				return false;
			}

		}
		return true;
	}

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new CompareBAMTagValues().instanceMain(args));
	}



}
