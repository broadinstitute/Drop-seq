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
import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.modularfileparser.DelimiterParser;
import org.broadinstitute.dropseqrna.utils.modularfileparser.ModularFileParser;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.RuntimeIOException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.nio.PicardHtsPath;

@CommandLineProgramProperties(summary = "Filters a BAM file based on a TAG and a file containing a list of values.  This is pretty similar to grepping with a file, but is faster and makes a proper BAM.", oneLineSummary = "Filters a BAM file based on a TAG and a file containing a list of values.", programGroup = DropSeq.class)
public class FilterBamByTag extends CommandLineProgram {

	private final Log log = Log.getInstance(FilterBamByTag.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze. This argument can accept wildcards, or a file with the suffix .bam_list that contains the locations of multiple BAM files", minElements = 1)
	public List<PicardHtsPath> INPUT;
	
	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output report")
	public File OUTPUT;
	
	@Argument(doc="A file containing a summary of the number of reads accepted and rejected.", optional=true)
	public File SUMMARY;
	
	@Argument(doc="Output a file containing the number of reads accepted and rejected for each tag.", optional=true)
	public File TAG_COUNTS_FILE;	

	@Argument(doc = "The BAM tag to filter on.")
	public String TAG;

	@Argument(doc = "A file with 1 column and 1 or more rows containing a barcode value per line.", optional = true)
	public File TAG_VALUES_FILE;

	@Argument(doc = "A tag value(s) for filtering reads.  Use instead of TAG_VALUES_FILE.", optional = true)
	public List<String> TAG_VALUE;

	@Argument(doc = "If having a tag value matches the values in the file, accept the read.  If set to false, reject the read.")
	public Boolean ACCEPT_TAG = true;
	
	@Argument(doc="Allow partial matching - if the tag value contains one of expected values, count as a match.  For example, if the allowed value is A and the tag is A,B, the read would match.", optional=true)
	public Boolean ALLOW_PARTIAL_MATCH=false;

	@Argument(doc = "In Paired Read Mode if the tag value is on either read the pair of reads is kept or discarded. This is slower when turned on because "
			+ "of the need to queryname sort the data, so only turn it on if you need it!")
	public Boolean PAIRED_MODE=false;

	@Argument(shortName="READ_MQ", doc = "Minimum mapping quality to include the read in the analysis.  Reads are not filtered on map quality by default.", optional=true)
	public Integer MINIMUM_MAPPING_QUALITY = null;

	@Argument(doc="If set to a a value < 1, the program will fail if fewer than this fraction of reads pass filters." +
			"  If set to a value >= 1, the program will fail if fewer than this many reads pass filters.", optional = true)
	public Double PASSING_READ_THRESHOLD;

	@Override
	protected int doWork() {
		if (TAG_VALUES_FILE == null && TAG_VALUE==null) {
			log.error("You must set either a file of tag values, or a single tag value.");
			return(1);
		}

		INPUT = FileListParsingUtils.expandPicardHtsPathList(INPUT);

		IOUtil.assertFileIsWritable(OUTPUT);
		if (this.SUMMARY!=null) IOUtil.assertFileIsWritable(this.SUMMARY);
		if (this.TAG_COUNTS_FILE!=null) IOUtil.assertFileIsWritable(this.TAG_COUNTS_FILE);
		Set<String> values;

		if (this.TAG_VALUES_FILE != null) {
			IOUtil.assertFileIsReadable(TAG_VALUES_FILE);
			values = readValues(this.TAG_VALUES_FILE);
		} else {
			values = new HashSet<>();
			if (this.TAG_VALUE!=null)
				values.addAll(this.TAG_VALUE);
		}

		final SamHeaderAndIterator headerAndIter =
				SamFileMergeUtil.mergeInputPaths(PicardHtsPath.toPaths(INPUT), false, SamReaderFactory.makeDefault());
		SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(headerAndIter.header, true, OUTPUT);

		if (!this.PAIRED_MODE)
			processUnpairedMode(headerAndIter, out, values, this.SUMMARY, this.TAG_COUNTS_FILE, MINIMUM_MAPPING_QUALITY, ALLOW_PARTIAL_MATCH);
		else
			processPairedMode(headerAndIter, out, values, this.SUMMARY, this.TAG_COUNTS_FILE, MINIMUM_MAPPING_QUALITY, ALLOW_PARTIAL_MATCH);

		return 0;
	}

	/**
	 * Just work through the reads one at a time.
	 * @param in
	 * @param out
	 */
	void processUnpairedMode (final SamHeaderAndIterator headerAndIter, final SAMFileWriter out, final Set<String> values, final File summaryFile, final File tagCountsFile, Integer mapQuality ,boolean allowPartial) {
		FilteredReadsMetric m = new FilteredReadsMetric();
		ProgressLogger progLog = new ProgressLogger(log);
		Iterator<SAMRecord> in = headerAndIter.iterator;
		while (in.hasNext()) {
			SAMRecord r = in.next();
			progLog.record(r);
			boolean filterFlag = filterRead(r, this.TAG, values, this.ACCEPT_TAG, mapQuality, allowPartial);
			if (!filterFlag) { 
				out.addAlignment(r);
				m.READS_ACCEPTED++;
				m.incrementTagAccepted(getTagValue(r.getAttribute(this.TAG)));				
			}
			else {
				m.READS_REJECTED++;
				m.incrementTagRejected(getTagValue(r.getAttribute(this.TAG)));
			}
		}
		writeSummary(summaryFile, m);
		writeTagCounts(tagCountsFile, m);
		CloserUtil.close(in);
		out.close();
		reportAndCheckFilterResults(m);
	}
	
	/**
	 * Write the summary output file of reads accepted/rejected.
	 * @param summaryOutputFile
	 * @param metrics
	 */
	private void writeSummary (File summaryOutputFile, FilteredReadsMetric metrics) {
		if (summaryOutputFile!=null) {
			MetricsFile<FilteredReadsMetric, Integer> outSummary = getMetricsFile();
			outSummary.addMetric(metrics);
			outSummary.write(summaryOutputFile);		
		}
	}
	
	private void writeTagCounts (File outFile, FilteredReadsMetric metrics) {
		if (outFile!=null) {
			MetricsFile<FilteredReadsMetric, String> outTagsCount = getMetricsFile();
			outTagsCount.addHistogram(metrics.getHistogramTagsAccepted());
			outTagsCount.addHistogram(metrics.getHistogramTagsRejected());
			outTagsCount.write(outFile);		
		}
	}

	/**
	 *
	 * @param in
	 * @param out
	 * @param values
	 */
	void processPairedMode (final SamHeaderAndIterator headerAndIter, final SAMFileWriter out, final Set<String> values, final File summaryFile, final File tagCountsFile, Integer mapQuality, boolean allowPartial) {
		ProgressLogger progLog = new ProgressLogger(log);
		FilteredReadsMetric m = new FilteredReadsMetric();
		
		PeekableIterator<SAMRecord> iter = new PeekableIterator<>(CustomBAMIterators.getQuerynameSortedRecords(headerAndIter));
		while (iter.hasNext()) {
			SAMRecord r1 = iter.next();
			progLog.record(r1);
			boolean filterFlag1 = filterRead(r1, this.TAG, values, this.ACCEPT_TAG, mapQuality, allowPartial);
			
			SAMRecord r2 = null;
			if (iter.hasNext()) r2 = iter.peek();
			// check for r2 being null in case the last read is unpaired.
			if (r2!=null && r1.getReadName().equals(r2.getReadName())) {
				// paired read found.
				progLog.record(r2);
				r2=iter.next();
				boolean filterFlag2 = filterRead(r2, this.TAG, values, this.ACCEPT_TAG, mapQuality, allowPartial);
				// if in accept tag mode, if either read shouldn't be filterd accept the pair
				// if in reject mode, if both reads shouldn't be filtered to accept the pair.
				if ((!filterFlag1 || !filterFlag2 & this.ACCEPT_TAG) || (!filterFlag1 && !filterFlag2 & !this.ACCEPT_TAG)) {
					out.addAlignment(r1);
					out.addAlignment(r2);
					m.READS_ACCEPTED+=2;
					m.incrementTagAccepted(getTagValue(r1.getAttribute(this.TAG)));
					m.incrementTagAccepted(getTagValue(r2.getAttribute(this.TAG)));
				} else { 
					m.READS_REJECTED+=2;
					m.incrementTagRejected(getTagValue(r1.getAttribute(this.TAG)));
					m.incrementTagRejected(getTagValue(r2.getAttribute(this.TAG)));
				}
			} else if (!filterFlag1) {
				out.addAlignment(r1);
				m.READS_ACCEPTED++;
				m.incrementTagAccepted(getTagValue(r1.getAttribute(this.TAG)));
			} else {
				m.READS_REJECTED++;
				m.incrementTagRejected(getTagValue(r1.getAttribute(this.TAG)));
			}
		}		
		writeSummary(summaryFile, m);
		writeTagCounts(tagCountsFile, m);
		out.close();
		reportAndCheckFilterResults(m);
	}

	private void reportAndCheckFilterResults(final FilteredReadsMetric m) {
		FilterProgramUtils.reportAndCheckFilterResults("reads", m.READS_ACCEPTED, m.READS_REJECTED,
				PASSING_READ_THRESHOLD, log);
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
			final boolean acceptFlag, final Integer mapQuality, boolean allowPartial) {

		// quickly filter on map quality if provided.
		if (mapQuality!=null && r.getMappingQuality() < mapQuality) return true;
		
		Object v = r.getAttribute(tag);
		// if the tag is not set, and you need a TAG set, then filter the read.
		if (v == null && acceptFlag)
			return true;
		// if the tag is not set, and you don't need a TAG set, then keep the
		// read.
		if (v == null && !acceptFlag)
			return false;

		final String vv =  getTagValue(v);

		// if there are no values to scan, it's a match. Start with that.
		boolean hasElement = true;
		
		// if there are values, check to see if this tag matches one.
		if (values != null && values.size() > 0)
			if (!allowPartial)
				hasElement = values.contains(vv);
			else {
				hasElement = values.stream().anyMatch((x -> vv.contains(x)));
			}

		if ((hasElement & acceptFlag)
				| (hasElement == false & acceptFlag == false))
			return false;
		if ((hasElement == false & acceptFlag)
				| (hasElement & acceptFlag == false))
			return true;

		return false;
	}
	
	private String getTagValue (Object v) {
		String vv = null;
		// short circuit if tag not set.
		if (v==null)
			return null;
		
		if (v instanceof Integer) {
			Integer o = (Integer) v;
			vv = Integer.toString(o);
		} else if (v instanceof String)
			vv = (String) v;
		else
			log.info("WHAT ELSE");
		return vv;
	}

	public static Set<String> readValues(File f) {
		try {
			f = f.getCanonicalFile();
		} catch (IOException e) {
			throw new RuntimeIOException("Exception reading " + f, e);
		}
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
		System.exit(new FilterBamByTag().instanceMain(args));
	}
}
