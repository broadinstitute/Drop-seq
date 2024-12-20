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
import java.io.PrintStream;
import java.util.*;

import htsjdk.samtools.*;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.SelectCellsByNumTranscripts;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;

import org.broadinstitute.dropseqrna.utils.alignmentcomparison.QueryNameJointIterator;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.*;
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

	@Argument(doc = "Output a version of INPUT_1 containing only the reads that are discordant with INPUT_2", optional = true)
	public File BAM_OUTPUT_1;

	@Argument(doc = "Output a version of INPUT_2 containing only the reads that are discordant with INPUT_1", optional = true)
	public File BAM_OUTPUT_2;

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

	@Argument(doc="If all tag values consist of the bases A,C,G,T (for example cell or molecular barcodes) then compress the internal representations to byte arrays to reduce memory usage." +
			"If your strings as of fixed length, use useFixedLengthBaseCompression instead, as it is even more memory efficient.  If values are encountered that can not be compressed the program will fail.", optional=true)
	public Boolean useBaseCompression=false;

	@Argument(doc="If all tag values consist of the bases A,C,G,T (for example cell or molecular barcodes) then compress the internal representations to byte arrays to reduce memory usage.  These strings must be of fixed length." +
			"If values are encountered that can not be compressed the program will fail.", optional=true)
	public Boolean useFixedLengthBaseCompression=false;

	@Argument(doc="Tags that do not have set values will be replaced with this value in the output. Default is '-'", optional=true)
	public String MISSING_TAG_VALUE="-";

	// a delimiter to use for the tag values if persisting as strings.
	private String TAG_DELIMITER=":";

	// Determine the mode
	private CompressionMode mode;

	// track the tag lengths for the fixed length compression mode.
	private int [] tagLengths;

	/**
	 * There are two ways to encode tag values for counting.
	 * 1: Keep the tag values as strings, use ":" as a delimiter between values as BAM tag values do not allow that string to be used.
	 * 2: Encode the tag values as byte arrays, where the byte array is the concatenation of the bytes of the tag values.  Values are encoded via Base64 encoding.
	 */
	ObjectCounter<String> tagValuesString = new ObjectCounter<>();
	ObjectCounter<ByteArrayWrapper> tagValuesByteArray = new ObjectCounter<>();

		@Override
	public int doWork() {
		// Expand the input file lists.
		this.INPUT_1 = FileListParsingUtils.expandPicardHtsPathList(INPUT_1);
		this.INPUT_2 = FileListParsingUtils.expandPicardHtsPathList(INPUT_2);

		// Determine the compression mode for output tags report.
		mode = CompressionMode.fromFlags(this.useBaseCompression, this.useFixedLengthBaseCompression);

		// Keep the headers to create output BAM files.
		final SamHeaderAndIterator headerAndIter_1 = SamFileMergeUtil.mergeInputPaths(
				PicardHtsPath.toPaths(this.INPUT_1), false, SamReaderFactory.makeDefault());

		final SamHeaderAndIterator headerAndIter_2 = SamFileMergeUtil.mergeInputPaths(
				PicardHtsPath.toPaths(this.INPUT_2), false, SamReaderFactory.makeDefault());

		// these values may be null if the output is not requested.
		SAMFileWriter writer_1 = createWriter(headerAndIter_1.header, this.BAM_OUTPUT_1);
		SAMFileWriter writer_2 = createWriter(headerAndIter_2.header, this.BAM_OUTPUT_2);

		PeekableIterator<List<SAMRecord>> iterator1 = getReadIterator (headerAndIter_1,this.READ_QUALITY, this.TAGS_1);
		PeekableIterator<List<SAMRecord>> iterator2 = getReadIterator (headerAndIter_2,this.READ_QUALITY, this.TAGS_2);

		log.info("Peakable Readers constructed");

		QueryNameJointIterator jointIterator = new QueryNameJointIterator(iterator1, iterator2);

		while (jointIterator.hasNext()) {
			QueryNameJointIterator.JointResult jr = jointIterator.next();
			// it's possible one data set is paired and, and the other is not.
			// given reads are filtered to have the same tag names, the tag values should be the same.
			// so it shouldn't matter which one is queried.
			List<SAMRecord> rList1 = jr.getOne();
			List<SAMRecord> rList2 = jr.getTwo();
			SAMRecord r1=rList1.getFirst();
			SAMRecord r2=rList2.getFirst();

			boolean tagsConcordant = compareTags(r1, r2, this.TAGS_1, this.TAGS_2);
			if (this.STRICT && !tagsConcordant) {
				log.error("Tag values differ for read: " + r1.getReadName());
				return 1;
			}
			writeDiscordantReads(writer_1, rList1, writer_2, rList2, tagsConcordant);
		}

		QueryNameJointIterator.QueryNameJointIteratorMetrics metrics= jointIterator.getMetrics();
		if (this.STRICT && metrics.hasDisjointReads()) {
			log.error("Input BAMs have different reads.  Exiting with error.");
			return 1;
		}

		if (this.READ_COUNT_OUTPUT!=null) {
			MetricsFile<QueryNameJointIterator.QueryNameJointIteratorMetrics, Integer> out = getMetricsFile();
			out.addMetric(jointIterator.getMetrics());
			out.write(READ_COUNT_OUTPUT);
		}

		writeTagValuesReport(TAG_VALUES_OUTPUT);

		log.info("DONE");
		return 0;
	}

	private void writeTagValuesReport (final File output) {
		if (this.TAG_VALUES_OUTPUT==null)
			return;

		PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(output));
		writeTagValuesHeader(out);
		switch (mode) {
			case NO_COMPRESSION -> {
				for (String key: tagValuesString.getKeys())
					out.println(String.join("\t", getTagValuesReportLine(key)));
			}
			case BASE_COMPRESSION, VARYING_LENGTH_COMPRESSION -> {
				for (ByteArrayWrapper key: tagValuesByteArray.getKeys())
					out.println(String.join("\t", getTagValuesReportLine(key)));
			}
		}
	}

	private void writeTagValuesHeader (PrintStream out) {
		List<String> header = new ArrayList<>();
		header.add("COUNT");
		for (int i=0; i<this.TAGS_1.size(); i++) {
			header.add("I1_"+this.TAGS_1.get(i));
			header.add("I2_"+this.TAGS_2.get(i));
		}
		out.println(String.join("\t", header));
	}

	private List<String> getTagValuesReportLine (String key) {
		List<String> result = new ArrayList<>();
		int count = tagValuesString.getCountForKey(key);
		String[] tagValues = key.split(this.TAG_DELIMITER);
		result.add(Integer.toString(count));
		result.addAll(Arrays.asList(tagValues));
		return (result);
	}

	private List<String> getTagValuesReportLine(ByteArrayWrapper key) {
    	List<String> tagValues;
    	switch (mode) {
        	case BASE_COMPRESSION -> tagValues = DNACompressor.decompressList(key.getData(), this.tagLengths);
        	case VARYING_LENGTH_COMPRESSION -> tagValues = DNACompressorVaryingLengths.decompressList(key.getData());
        	default -> throw new IllegalStateException("Unexpected compression mode: " + mode);
    	}
    	return tagValues;
	}

	private void writeDiscordantReads (final SAMFileWriter writer1, List<SAMRecord> r1, final SAMFileWriter writer2, List<SAMRecord> r2, boolean tagsConcordant) {
		if (tagsConcordant)
			return;
		if (writer1!=null)
			r1.forEach(writer1::addAlignment);
		if (writer2!=null)
			r2.forEach(writer2::addAlignment);
	}

	private SAMFileWriter createWriter(SAMFileHeader header, File output) {
		return output != null ? new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, output) : null;
	}

	private boolean compareTags (final SAMRecord r1, final SAMRecord r2, final List<String> tags1, final List<String> tags2) {
		boolean sameRead = r1.getReadName().equals(r2.getReadName());
		// this should never happen.
		if (!sameRead)
			log.error("Read names differ at iteration.  R1: "+ r1.toString(), "R2: ", r2.toString());

		boolean tagValuesConcordant = true;
		List<String> tagValuesList = new ArrayList<>();

		for (int i=0; i<tags1.size(); i++) {
			String tag1 = tags1.get(i);
			String tag2 = tags2.get(i);
			String v1 = getTagValue(tag1, r1);
			String v2 = getTagValue(tag2, r2);

			if (!v1.equals(v2))
				tagValuesConcordant=false;
			tagValuesList.add(v1);
			tagValuesList.add(v2);
		}

		saveTagValues(tagValuesList, tagValuesConcordant);

		return tagValuesConcordant;
	}

	/**
	 * Gets the value for the tag for this read.  If the tag is not set, return the missing tag value.
	 * @return The value of the tag, or the missing tag value if the tag is not set.
	 */
	private String getTagValue (String tag, SAMRecord r) {
		String result = r.getStringAttribute(tag);
		if (result==null)
			result=this.MISSING_TAG_VALUE;
		return result;
	}

	private void saveTagValues (List<String> tagValuesList, boolean tagValuesConcordant) {
		if (DISCORDANT_READS_ONLY && tagValuesConcordant)
			return;

		switch (mode) {
			case NO_COMPRESSION -> {
				String tagValues = String.join(TAG_DELIMITER, tagValuesList);
				tagValuesString.increment(tagValues);
			}
			case BASE_COMPRESSION -> {
				byte[] tagValuesBytes = DNACompressor.compressList(tagValuesList);
				// special case, where tag values must have consistent lengths.
				if (!validateTagLengths(tagValuesList))
					throw new IllegalArgumentException("Tag values have inconsistent lengths, but fixed length compression is enabled. Tag values: " + tagValuesList);
				tagValuesByteArray.increment(new ByteArrayWrapper(tagValuesBytes));
			}
			case VARYING_LENGTH_COMPRESSION -> {
				byte[] tagValuesBytes = DNACompressorVaryingLengths.compressList(tagValuesList);
				tagValuesByteArray.increment(new ByteArrayWrapper(tagValuesBytes));
			}
		}
	}

	private boolean validateTagLengths(List<String> tagValuesList) {
		// Compute the lengths of the current tag values
		int[] lengths = tagValuesList.stream()
				.mapToInt(String::length)
				.toArray();

		// Initialize tagLengths if not already set
		if (this.tagLengths == null) {
			this.tagLengths = lengths;
			return true; // Always valid on first initialization
		}

		// Validate the current lengths against the stored tagLengths
		return Arrays.equals(lengths, this.tagLengths);
	}

	private PeekableIterator<List<SAMRecord>> getReadIterator (SamHeaderAndIterator headerAndIter, final Integer readQuality, List<String> requiredTags) {
		// sort the data by query name.
		// Iterator<SAMRecord> iter = getQueryNameSortedData(headerAndIter);
		Iterator<SAMRecord> iter = headerAndIter.iterator;

		// Filter out reads below a map quality threshold.  This removes non-primary reads.
		iter = new MapQualityFilteredIterator(iter, readQuality, true).iterator();

		// drop reads that do not have the required tags set
		iter = new MissingTagFilteringIterator(iter, requiredTags.toArray(new String[0]));

		// Queryname sort data if it is not already sorted by queryname.  This will spill to disk.
		// if the data was already in queryname order, it can be streamed and filtered directly.
		if (!headerAndIter.header.getSortOrder().equals(SAMFileHeader.SortOrder.queryname)) {
			log.info("Input SAM/BAM not in queryname order, sorting...");
			final ProgressLogger progressLogger = new ProgressLogger(log, 1000000, "Sorting reads in query name order");
			iter = SamRecordSortingIteratorFactory.create(headerAndIter.header, headerAndIter.iterator, READ_NAME_COMPARATOR, progressLogger);
		}

		final GroupingIterator<SAMRecord> groupingIterator = new GroupingIterator<>(iter, READ_NAME_COMPARATOR);
		PeekableIterator<List<SAMRecord>> peekableIterator = new PeekableIterator<>(groupingIterator);
		return peekableIterator;
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

	public enum CompressionMode {
		NO_COMPRESSION(false, false),
		BASE_COMPRESSION(true, false),
		VARYING_LENGTH_COMPRESSION(false, true);

		private final boolean useBaseCompression;
		private final boolean useFixedLengthBaseCompression;

		// Constructor
		CompressionMode(boolean useBaseCompression, boolean useFixedLengthBaseCompression) {
			this.useBaseCompression = useBaseCompression;
			this.useFixedLengthBaseCompression = useFixedLengthBaseCompression;
		}

		// Method to determine the mode based on flags
		public static CompressionMode fromFlags(boolean useBaseCompression, boolean useFixedLengthBaseCompression) {
			for (CompressionMode mode : CompressionMode.values()) {
				if (mode.useBaseCompression == useBaseCompression && mode.useFixedLengthBaseCompression == useFixedLengthBaseCompression) {
					return mode;
				}
			}
			throw new IllegalArgumentException("Invalid flag combination");
		}
	}


	protected String[] customCommandLineValidation() {

		final ArrayList<String> list = new ArrayList<>(1);


		if (this.TAGS_1.size()!=this.TAGS_2.size())
			list.add("TAGS_1 and TAGS_2 must be the same length.");

		if (BAM_OUTPUT_1 != null)
			IOUtil.assertFileIsReadable(this.BAM_OUTPUT_1);

		if (BAM_OUTPUT_2 != null)
			IOUtil.assertFileIsReadable(this.BAM_OUTPUT_2);

		if (READ_COUNT_OUTPUT != null)
			IOUtil.assertFileIsReadable(this.READ_COUNT_OUTPUT);

		if (TAG_VALUES_OUTPUT != null)
			IOUtil.assertFileIsReadable(this.TAG_VALUES_OUTPUT);

		if (useBaseCompression & useFixedLengthBaseCompression) {
			list.add("Overriding the useBaseCompression flag with useFixedLengthBaseCompression flag.");
			this.useBaseCompression=false;
		}



		return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);
	}





	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new CompareBAMTagValues().instanceMain(args));
	}



}
