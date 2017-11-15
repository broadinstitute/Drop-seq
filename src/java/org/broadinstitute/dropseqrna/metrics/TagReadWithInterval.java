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

import java.io.File;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.broadinstitute.dropseqrna.cmdline.DropSeq;

import org.broadinstitute.dropseqrna.utils.SamHeaderUtil;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(usage = "", usageShort = "", programGroup = DropSeq.class)
public class TagReadWithInterval extends CommandLineProgram {
	private final Log log = Log.getInstance(TagReadWithInterval.class);
	
	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.  Must be coordinate sorted.")
	public File INPUT;

	@Option(doc = "The list of Loci to gather start/end read positons for.  This file is in Interval format - tab seperated with fields: chr start end strand name")
	public File LOCI;

	@Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output BAM, written with new interval tag")
	public File OUTPUT;
	
	@Option(doc = "The tag name to use.  Defaults to ZI.  If a read previously had a tag and now does not, the tag is removed.", optional=true)
	public String TAG="ZI";

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new TagReadWithInterval().instanceMain(args));
	}
	
	@Override
	protected int doWork() {
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		SamReader inputSam = SamReaderFactory.makeDefault().open(INPUT);
		SAMFileHeader header = inputSam.getFileHeader();
        header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
		SamHeaderUtil.addPgRecord(header, this);

		SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT);
        
		IntervalList loci = IntervalList.fromFile(this.LOCI);
		OverlapDetector<Interval> od =getOverlapDetector (loci);
		ProgressLogger processLogger = new ProgressLogger(log);
		
		for (SAMRecord record: inputSam) {
			processLogger.record(record);
			
			if (record.getReadUnmappedFlag()==false) {
				// use alignment blocks instead of start/end to properly deal with split reads mapped over exon/exon boundaries.
				record=tagRead(record, od);				
			} else {
				record.setAttribute(this.TAG, null);
			}
			writer.addAlignment(record);
			
		}
		
		CloserUtil.close(inputSam);
		writer.close();
		
		return(0);
		
		
	}
	
	private SAMRecord tagRead(SAMRecord r, OverlapDetector<Interval> od) {
		Set<Interval>intervals = new HashSet<Interval>();
		List<AlignmentBlock> blocks = r.getAlignmentBlocks();
		for (AlignmentBlock b: blocks) {
			int refStart =b.getReferenceStart();
			int refEnd = refStart+b.getLength()-1;
			Interval v = new Interval(r.getReferenceName(), refStart, refEnd);
			Collection<Interval> blockResult = od.getOverlaps(v);
			intervals.addAll(blockResult);
		}
		String tagName = getIntervalName(intervals);
		if (tagName!=null) {
			r.setAttribute(this.TAG, tagName);
		} else {
			r.setAttribute(this.TAG, null);
		}
		return (r);
		
	}
	
	/**
	 * 
	 * @param intervals
	 * @return Returns a comma seperated list of interval names, or null if the list is empty.
	 */
	private String getIntervalName (Collection<Interval> intervals) {
		
		if (intervals.isEmpty()) return (null);
		
		StringBuilder result = new StringBuilder();
		Iterator<Interval> iter = intervals.iterator();
		result.append(iter.next().getName());
		
		
		while (iter.hasNext()) {
			result.append(",");
			result.append(iter.next().getName());
		}
		
		return (result.toString());
	}
	
	
	
	/**
	 * Each interval has a corresponding ReadDepthMetric.
	 * @param loci
	 * @return
	 */
	OverlapDetector<Interval> getOverlapDetector (IntervalList loci) {
		OverlapDetector<Interval> od = new OverlapDetector<Interval>(0, 0);
		
		for (Interval i: loci.getIntervals()) {
			od.addLhs(i, i);
		}
		return (od);
	}
	
}
