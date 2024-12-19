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
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.SelectCellsByNumTranscripts;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import org.broadinstitute.dropseqrna.utils.alignmentcomparison.QueryNameJointIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.MapQualityFilteredIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamRecordSortingIteratorFactory;
import picard.cmdline.CommandLineProgram;
import picard.nio.PicardHtsPath;

@CommandLineProgramProperties(summary = "Tests that two BAMs have the same TAG values per read.  "
		+ "We assume the two BAMs tested have the same set of reads in the same order, and test the read names "
		+ "for each iteration match, and the values for the TAGs listed are equal.",
oneLineSummary = "Tests that two BAMs have the same TAG values per read.",
programGroup = DropSeq.class)

public class CompareBAMTagValues extends CommandLineProgram {

	private final Log log = Log.getInstance(CompareBAMTagValues.class);

	@Argument(doc = "The input SAM or BAM file to analyze.", optional=false)
	public List<PicardHtsPath> INPUT_1;

	@Argument(doc = "The input SAM or BAM file to analyze.", optional=false)
	public List<PicardHtsPath> INPUT_2;

	@Argument(doc="A list of tags to test in INPUT_1.  These tags are paired with TAGS_2 at the same index.  For each index, the values of TAGS_1[i] in INPUT_1 are" +
			"compared to the values of TAGS_2[i] in INPUT_2 for equality.", optional=false)
	public List<String> TAGS_1;

	@Argument(doc="A list of tags to test in INPUT_2.  Values of these tags will be compared across BAMs for equality.", optional=false)
	public List<String> TAGS_2;

	@Argument(doc="Output file of the tag values that differ between the two BAMs.  The program records each unique tuple of tag values with the number of reads that contain those values." +
			"The column names contain the input identifier I1 or I2, followed by the TAGS_1 and TAGS_2 names for each index.", optional=true)
	public File TAG_VALUES_OUTPUT;

	@Argument(doc="If true, only discordant read tag values are tracked and emitted.  For tags with many possible values, this option can be more memory efficient.  " +
			"If false, all read tag values are emitted.", optional=true)
	public Boolean DISCORDANT_READS_ONLY=false;

	@Argument(doc="A summary of the number of reads that were observed in each BAM.  Records the number of reads only in INPUT_1, only in INPUT_2, and in both INPUT_1 and INPUT_2.", optional=true)
	public File READ_COUNT_OUTPUT;

	@Argument(doc="If true, the program will exit with an error if any read is encountered that is not in both INPUTS or has TAG values that disagree")
	public Boolean STRICT=true;

	@Argument(doc="Restrict analysis to reads with a mapping quality greater than or equal to this value.  Reads with a mapping quality less than this value are ignored.  Non-primary reads are always ignored.", optional=true)
	public int READ_QUALITY=0;

	@Argument(doc="If all tag values consist of the bases A,C,G,T (for example cell or molecular barcodes) then compress the internal representations to byte arrays to reduce memory usage.  These strings must be of fixed length." +
			"If values are encountered that can not be compressed the program will fail.", optional=true)
	public Boolean useBaseCompression=false;

	@Override
	public int doWork() {
		this.INPUT_1 = FileListParsingUtils.expandPicardHtsPathList(INPUT_1);
		this.INPUT_2 = FileListParsingUtils.expandPicardHtsPathList(INPUT_2);

		// Keep the headers to create output BAM files.
		final SamHeaderAndIterator headerAndIter_1 = SamFileMergeUtil.mergeInputPaths(
				PicardHtsPath.toPaths(this.INPUT_1), false, SamReaderFactory.makeDefault());

		final SamHeaderAndIterator headerAndIter_2 = SamFileMergeUtil.mergeInputPaths(
				PicardHtsPath.toPaths(this.INPUT_1), false, SamReaderFactory.makeDefault());


		PeekableIterator<List<SAMRecord>> iterator1 = getReadIterator (headerAndIter_1,READ_QUALITY);
		PeekableIterator<List<SAMRecord>> iterator2 = getReadIterator (headerAndIter_2,READ_QUALITY);

		QueryNameJointIterator jointIterator = new QueryNameJointIterator(iterator1, iterator2);

		while (jointIterator.hasNext()) {
			QueryNameJointIterator.JointResult jr = jointIterator.next();
			List<SAMRecord> r1 = jr.getOne();
			List<SAMRecord> r2 = jr.getTwo();


		}


		MetricsFile<QueryNameJointIterator.QueryNameJointIteratorMetrics, Integer> out = getMetricsFile();
		out.addMetric(jointIterator.getMetrics());
		out.write(READ_COUNT_OUTPUT);

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

	private PeekableIterator<List<SAMRecord>> getReadIterator (SamHeaderAndIterator headerAndIter, final Integer readQuality) {
		// sort the data by query name.
		Iterator<SAMRecord> iter = getQueryNameSortedData(headerAndIter);

		// optionally, filter out reads below a map quality threshold.
		if (readQuality!=null) iter = new MapQualityFilteredIterator(iter, readQuality, true).iterator();

		//TODO: Do I need to group reads?  Maybe if I want to be tollerent of paired read data.
		final GroupingIterator<SAMRecord> groupingIterator = new GroupingIterator<>(iter, READ_NAME_COMPARATOR);
		PeekableIterator<List<SAMRecord>> peekableIterator = new PeekableIterator<>(groupingIterator);
		return peekableIterator;
	}

	/**
	 * Sort the data by query name if it is not already sorted by query name.
	 * @param headerAndIter The input data.
	 * @return An iterator sorted by queryname.
	 */
	private Iterator<SAMRecord> getQueryNameSortedData (final SamHeaderAndIterator headerAndIter) {
		if (headerAndIter.header.getSortOrder().equals(SAMFileHeader.SortOrder.queryname))
			return headerAndIter.iterator;
		log.info("Input SAM/BAM not in queryname order, sorting...");
		final ProgressLogger progressLogger = new ProgressLogger(log, 1000000, "Sorting reads in query name order");
		final CloseableIterator<SAMRecord> result = SamRecordSortingIteratorFactory.create(headerAndIter.header, headerAndIter.iterator, READ_NAME_COMPARATOR, progressLogger);
		log.info("Sorting finished.");
		return result;
	}

	/**
	 * A comparator that sorts reads by query name.
	 * This compares more than just the query name.
	 * @return A metrics file.
	 */
	static final Comparator<SAMRecord> READ_NAME_COMPARATOR =  new Comparator<SAMRecord>() {
		private final SAMRecordQueryNameComparator comp = new SAMRecordQueryNameComparator();
		@Override
		public int compare(final SAMRecord s1, final SAMRecord s2) {
			return comp.fileOrderCompare(s1, s2);
		}
	};


	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new CompareBAMTagValues().instanceMain(args));
	}



}
