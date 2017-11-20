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

import java.io.File;
import java.util.HashSet;
import java.util.Set;

import htsjdk.samtools.util.*;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.modularfileparser.DelimiterParser;
import org.broadinstitute.dropseqrna.utils.modularfileparser.ModularFileParser;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(usage = "Filters a BAM file based on a TAG and a file containing a list of values.  This is pretty similar to grepping with a file, but is faster and makes a proper BAM.", usageShort = "Filters a BAM file based on a TAG and a file containing a list of values.", programGroup = DropSeq.class)
public class FilterBAMByTag extends CommandLineProgram {

	private final Log log = Log.getInstance(FilterBAMByTag.class);

	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.")
	public File INPUT;

	@Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output report")
	public File OUTPUT;

	@Option(doc = "The BAM tag to filter on.")
	public String TAG;

	@Option(doc = "A file with 1 column and 1 or more rows containing a barcode value per line.", optional = true)
	public File TAG_VALUES_FILE;

	@Option(doc = "A single value for filtering reads.  Use instead of TAG_VALUES_FILE.", optional = true)
	public String TAG_VALUE;

	@Option(doc = "If having a tag value matches the values in the file, accept the read.  If set to false, reject the read.")
	public Boolean ACCEPT_TAG = true;

	@Option(doc = "In Paired Read Mode if the tag value is on either read the pair of reads is kept or discarded. This is slower when turned on because "
			+ "of the need to queryname sort the data, so only turn it on if you need it!")
	public Boolean PAIRED_MODE=false;

	@Override
	protected int doWork() {
		if (TAG_VALUES_FILE == null && TAG == null) {
			log.error("You must set either a file of tag values, or a single tag value.");
			System.exit(1);
		}

		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		Set<String> values;

		if (this.TAG_VALUES_FILE != null) {
			IOUtil.assertFileIsReadable(TAG_VALUES_FILE);
			values = readValues(this.TAG_VALUES_FILE);
		} else {
			values = new HashSet<>();
			if (this.TAG_VALUE!=null)
				values.add(this.TAG_VALUE);
		}

		SamReader in = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE).open(INPUT);
		SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(
				in.getFileHeader(), true, OUTPUT);

		if (!this.PAIRED_MODE)
			processUnpairedMode(in, out, values);
		else
			processPairedMode(in, out, values);

		return 0;
	}

	/**
	 * Just work through the reads one at a time.
	 * @param in
	 * @param out
	 */
	void processUnpairedMode (final SamReader in, final SAMFileWriter out, final Set<String> values) {
		ProgressLogger progLog = new ProgressLogger(log);
		for (final SAMRecord r : in) {
			progLog.record(r);
			boolean filterFlag = filterRead(r, this.TAG, values, this.ACCEPT_TAG);
			if (!filterFlag)
				out.addAlignment(r);
		}
		CloserUtil.close(in);
		out.close();
	}

	/**
	 *
	 * @param in
	 * @param out
	 * @param values
	 */
	void processPairedMode (final SamReader in, final SAMFileWriter out, final Set<String> values) {
		ProgressLogger progLog = new ProgressLogger(log);
		PeekableIterator<SAMRecord> iter = new PeekableIterator<>(CustomBAMIterators.getQuerynameSortedRecords(in));
		while (iter.hasNext()) {
			SAMRecord r1 = iter.next();
			progLog.record(r1);
			boolean filterFlag1 = filterRead(r1, this.TAG, values, this.ACCEPT_TAG);

			SAMRecord r2 = null;
			if (iter.hasNext()) r2 = iter.peek();
			// check for r2 being null in case the last read is unpaired.
			if (r2!=null && r1.getReadName().equals(r2.getReadName())) {
				// paired read found.
				progLog.record(r2);
				r2=iter.next();
				boolean filterFlag2 = filterRead(r2, this.TAG, values, this.ACCEPT_TAG);
				// if in accept tag mode, if either read shouldn't be filterd accept the pair
				// if in reject mode, if both reads shouldn't be filtered to accept the pair.
				if ((!filterFlag1 || !filterFlag2 & this.ACCEPT_TAG) || (!filterFlag1 && !filterFlag2 & !this.ACCEPT_TAG)) {
					out.addAlignment(r1);
					out.addAlignment(r2);
				}
			} else if (!filterFlag1)
				out.addAlignment(r1);
		}
		CloserUtil.close(in);
		out.close();
	}

	boolean retainByReadNumber (final SAMRecord r, final int desiredReadNumber) {
		boolean readPaired = r.getReadPairedFlag();
		// if the read isn't paired, then there's no read number, so keep it.
		if (!readPaired)
			return true;
		// read is paired, is it the first of pair?
		boolean firstRead = r.getFirstOfPairFlag();
		// if you're looking for the first read and this isn't, or looking for the 2nd read and this isn't, then go to the next read.
		if ((firstRead && desiredReadNumber!=1) || (!firstRead && desiredReadNumber==1)) return (false);
		return true;
	}

	boolean filterRead(final SAMRecord r, final String tag, final Set<String> values,
			final boolean acceptFlag) {
		Object v = r.getAttribute(tag);
		// if the tag is not set, and you need a TAG set, then filter the read.
		if (v == null && acceptFlag)
			return true;
		// if the tag is not set, and you don't need a TAG set, then keep the
		// read.
		if (v == null && !acceptFlag)
			return false;

		String vv = null;

		if (v instanceof Integer) {
			Integer o = (Integer) v;
			vv = Integer.toString(o);
		} else if (v instanceof String)
			vv = (String) v;
		else
			log.info("WHAT ELSE");

		// if there are no values to scan, it's a match. Start with that.
		boolean hasElement = true;

		// if there are values, check to see if this tag matches one.
		if (values != null && values.size() > 0)
			hasElement = values.contains(vv);

		if ((hasElement & acceptFlag)
				| (hasElement == false & acceptFlag == false))
			return false;
		if ((hasElement == false & acceptFlag)
				| (hasElement & acceptFlag == false))
			return true;

		return false;
	}

	private Set<String> readValues(final File f) {
		Set<String> result = new HashSet<>();
		ModularFileParser p = new ModularFileParser(new DelimiterParser(","),
				f, 0);
		String[] line = null;
		while ((line = p.readNextLine()) != null)
			result.add(line[0]);
		p.close();
		return result;
	}

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new FilterBAMByTag().instanceMain(args));
	}
}
