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
package org.broadinstitute.dropseqrna.metrics;

import java.io.File;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.SamHeaderUtil;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.samtools.util.ProgressLogger;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(summary = "", oneLineSummary = "", programGroup = DropSeq.class)
public class TagReadWithInterval extends CommandLineProgram {
	private final Log log = Log.getInstance(TagReadWithInterval.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.  Must be coordinate sorted.")
	public File INPUT;

	@Argument(doc = "The list of intervals to tag onto reads in the BAM.  This file is in Interval format - tab seperated with fields: chr start end strand name")
	public File INTERVALS;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output BAM, written with new interval tag")
	public File OUTPUT;

	@Argument(doc = "The tag name to use.  Defaults to ZI.  If a read previously had a tag and now does not, the tag is removed.", optional=true)
	public String TAG="ZI";

	@Override
	protected int doWork() {
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		SamReader inputSam = SamReaderFactory.makeDefault().open(INPUT);
		SAMFileHeader header = inputSam.getFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
		SamHeaderUtil.addPgRecord(header, this);

		SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT);

		IntervalList loci = IntervalList.fromFile(this.INTERVALS);
		OverlapDetector<Interval> od =getOverlapDetector (loci);
		ProgressLogger processLogger = new ProgressLogger(log);

		for (SAMRecord record: inputSam) {
			processLogger.record(record);

			if (record.getReadUnmappedFlag()==false)
				// use alignment blocks instead of start/end to properly deal with split reads mapped over exon/exon boundaries.
				record=tagRead(record, od);
			else
				record.setAttribute(this.TAG, null);
			writer.addAlignment(record);

		}

		CloserUtil.close(inputSam);
		writer.close();

		return(0);


	}

	private SAMRecord tagRead(final SAMRecord r, final OverlapDetector<Interval> od) {
		Set<Interval>intervals = new HashSet<>();
		List<AlignmentBlock> blocks = r.getAlignmentBlocks();
		for (AlignmentBlock b: blocks) {
			int refStart =b.getReferenceStart();
			int refEnd = refStart+b.getLength()-1;
			Interval v = new Interval(r.getReferenceName(), refStart, refEnd);
			Collection<Interval> blockResult = od.getOverlaps(v);
			intervals.addAll(blockResult);
		}
		String tagName = getIntervalName(intervals);
		if (tagName!=null)
			r.setAttribute(this.TAG, tagName);
		else
			r.setAttribute(this.TAG, null);
		return (r);

	}

	/**
	 *
	 * @param intervals
	 * @return Returns a comma separated list of interval names, or null if the list is empty.
	 */
	String getIntervalName (final Collection<Interval> intervals) {
		if (intervals.isEmpty()) return (null);
		List<String> names = intervals.stream().map(x-> x.getName()).collect(Collectors.toList());
		String result = StringUtils.join(names, ",");
		return result;

	}


	/**
	 * Each interval has a corresponding ReadDepthMetric.
	 * @param loci
	 * @return
	 */
	OverlapDetector<Interval> getOverlapDetector (final IntervalList loci) {
		OverlapDetector<Interval> od = new OverlapDetector<>(0, 0);

		for (Interval i: loci.getIntervals())
			od.addLhs(i, i);
		return (od);
	}

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new TagReadWithInterval().instanceMain(args));
	}


}
