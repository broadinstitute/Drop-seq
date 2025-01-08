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

import java.util.regex.Pattern;
import htsjdk.samtools.*;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;

import org.broadinstitute.dropseqrna.utils.alignmentcomparison.QueryNameJointIterator;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.*;
import picard.cmdline.CommandLineProgram;
import picard.nio.PicardHtsPath;

import java.util.stream.IntStream;
import java.util.stream.Collectors;

@CommandLineProgramProperties(summary = "Tests that two BAMs have the same TAG values per read.  "
		+ "This program only considers primary reads, and can optional filter to reads with a minimum mapping quality."
		+ "Reads are matched by query name string and first or second of the pair (See Picard's SAMRecordQueryNameComparator)."
		+ "Tags on reads are then compared for equality.  If tag values differ or tag values are missing, the read is considered discordant.",
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

	@Argument(doc="Some programs (STARsolo) set missing tag values to '-'.  Tags with this value are treated as unset for comparisons.  Set to null to disable this filter.", optional=true)
	public String TAG_MISSING_VALUE = "-";

	@Argument(doc="Some programs (CellRanger) have suffixes on tags like the cell barcode. This option will trim that suffix from all tag values", optional=true)
	public String REMOVE_TAG_SUFFIX="-1";

	@Argument(doc="Output file of the tag values that differ between the two BAMs.  The program records each unique tuple of tag values with the number of reads that contain those values." +
			"The column names contain the input identifier I1 or I2, followed by the TAGS_1 and TAGS_2 names for each index.", optional=true)
	public File TAG_VALUES_OUTPUT;

	@Argument(doc="If true, only discordant read tag values are tracked and emitted in TAG_VALUES_OUTPUT.  For tags with many possible values, this option can be more memory efficient.  " +
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
	private Boolean useBaseCompression=false;

	@Argument(doc="If all tag values consist of the bases A,C,G,T (for example cell or molecular barcodes) then compress the internal representations to byte arrays to reduce memory usage.  These strings must be of fixed length." +
			"If values are encountered that can not be compressed the program will fail.", optional=true)
	public Boolean useFixedLengthBaseCompression=false;

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
		if (TAG_VALUES_OUTPUT!=null)
			log.info("Using tag value compression mode: " + mode);

		// Keep the headers to create output BAM files.
		final SamHeaderAndIterator headerAndIter_1 = SamFileMergeUtil.mergeInputPaths(
				PicardHtsPath.toPaths(this.INPUT_1), false, SamReaderFactory.makeDefault());

		final SamHeaderAndIterator headerAndIter_2 = SamFileMergeUtil.mergeInputPaths(
				PicardHtsPath.toPaths(this.INPUT_2), false, SamReaderFactory.makeDefault());

		// these values may be null if the output is not requested.
		SAMFileWriter writer_1 = createWriter(headerAndIter_1.header, this.BAM_OUTPUT_1);
		SAMFileWriter writer_2 = createWriter(headerAndIter_2.header, this.BAM_OUTPUT_2);

		log.info("Constructing BAM reader for input" + this.INPUT_1.getFirst().toString());
		PeekableIterator<List<SAMRecord>> iterator1 = getReadIterator (headerAndIter_1,this.READ_QUALITY, this.TAGS_1);

		log.info("Constructing BAM reader for input" + this.INPUT_2.getFirst().toString());
		PeekableIterator<List<SAMRecord>> iterator2 = getReadIterator (headerAndIter_2,this.READ_QUALITY, this.TAGS_2);

		QueryNameJointIterator jointIterator = new QueryNameJointIterator(iterator1, iterator2);

		while (jointIterator.hasNext()) {
			QueryNameJointIterator.JointResult jr = jointIterator.next();
			// it's possible one data set is paired and, and the other is not.
			// if one read is paired and the other is not, then if either read of the pair is concordant, the pair is concordant.
			// if both reads are paired then compare both reads and mark as discordant if either read has discordant tags.
			List<SAMRecord> rList1 = jr.getOne();
			List<SAMRecord> rList2 = jr.getTwo();
			String readName = rList1.getFirst().getReadName();

			boolean tagsConcordant = compareTagsReadGroup(rList1, rList2, this.TAGS_1, this.TAGS_2);
			if (this.STRICT && !tagsConcordant) {
				log.error("Tag values differ for read: " + readName);
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

		// close the writers.
		if (writer_1!=null)
			writer_1.close();
		if (writer_2!=null)
			writer_2.close();

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
				// Sort keys alphanumerically
				List<String> sortedKeys = new ArrayList<>(tagValuesString.getKeys());
				sortedKeys.sort(CustomKeyComparator.INSTANCE); // Apply custom sorting by components
				// Replace "null" with "NA" after sorting
				// sortedKeys.replaceAll(key -> key == null ? "NA" : key.replaceAll("\\bnull\\b", "NA"));

				for (String key : sortedKeys) {
					// Get the report line and apply the null -> "NA" transformation
					List<String> reportLine = replaceNullsWithNA(getTagValuesReportLine(key));
					out.println(String.join("\t", reportLine));
				}
			}
			case FIXED_LENGTH_COMPRESSION, VARYING_LENGTH_COMPRESSION -> {
				// Sort keys by compressed output
				List<ByteArrayWrapper> sortedKeys = new ArrayList<>(tagValuesByteArray.getKeys());
				sortedKeys.sort(Comparator.comparing(
						ByteArrayWrapper::getData, // Extract raw byte arrays from ByteArrayWrapper
						DNACompressor.getDecompressedComparator(this.tagLengths) // Use comparator from DNACompressor
				));

				for (ByteArrayWrapper key : sortedKeys) {
					// Get the report line and apply the null -> "NA" transformation
					List<String> reportLine = replaceNullsWithNA(getTagValuesReportLine(key));
					out.println(String.join("\t", reportLine));
				}
			}
		}

	}

	// Reusable method to replace null values with "NA" in a list
	private static List<String> replaceNullsWithNA(List<String> list) {
		if (list != null) {
			list.replaceAll(value -> {
				if (value == null) {
					return "NA"; // Replace completely null values
				}
				// Replace "null" substrings in concatenated strings with "NA"
				return value.replaceAll("\\bnull\\b", "NA");
			});
		}
		return list;
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
        	case FIXED_LENGTH_COMPRESSION -> tagValues = DNACompressor.decompressList(key.getData(), this.tagLengths);
        	case VARYING_LENGTH_COMPRESSION -> tagValues = DNACompressorVaryingLengths.decompressList(key.getData());
        	default -> throw new IllegalStateException("Unexpected compression mode: " + mode);
    	}
		List<String> result = new ArrayList<>();
		int count = tagValuesByteArray.getCountForKey(key);
		result.add(Integer.toString(count));
		result.addAll(tagValues);
		return (result);
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
		// if I don't clone the header, the side effect breaks the main program iterator.
		SAMFileHeader outputHeader = header.clone();
		outputHeader.setSortOrder(SAMFileHeader.SortOrder.queryname);
		return output != null ? new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, output) : null;
	}

	/**
	 * Compares the tags of reads between two lists to determine if they are concordant.
	 *
	 * Logic:
	 * - If both lists have a single read, compare the tags directly using {@code compareTags}.
	 * - If one list has a single read, compare it against all reads in the other (multi-read) list:
	 *     - Filter the multi-read list to remove reads missing the required tags.
	 *     - If no matching reads are found, the first read from the filtered list is selected.
	 *     - Determine if any read in the multi-read list is concordant with the single read.
	 * - If both lists have multiple reads:
	 *     - Verify that both lists have the same size; if not, an exception is thrown.
	 *     - Sort both lists by first/second of pair using {@code PAIRED_READ_ORDER_COMPARATOR}.
	 *     - Compare the tags of reads pairwise. All reads must match for the result to be true.
	 *
	 * Behavior:
	 * - Returns {@code true} if all relevant reads are concordant, {@code false} otherwise.
	 * - Throws an {@link IllegalStateException} if both lists have multiple reads but different sizes.
	 *
	 * @param r1List The first list of reads. This can be a single-read or multi-read list.
	 * @param r2List The second list of reads. This can also be a single-read or multi-read list.
	 * @param tags1 The tags to compare in the reads from {@code r1List}.
	 * @param tags2 The tags to compare in the reads from {@code r2List}.
	 * @return {@code true} if the tags are concordant across the two lists; {@code false} otherwise.
	 */
	private boolean compareTagsReadGroup(final List<SAMRecord> r1List,  final List<SAMRecord> r2List,  final List<String> tags1,  final List<String> tags2) {
		// Handle cases with one list empty (optional, depending on your requirements)
		if (r1List.isEmpty() || r2List.isEmpty()) {
			return false;
		}

		// Both lists have a single read
		if (r1List.size() == 1 && r2List.size() == 1) {
			return compareTags(r1List.getFirst(), r2List.getFirst(), tags1, tags2, true);
		}

		// One list has a single read, compare against all reads in the other list
		if (r1List.size() == 1 || r2List.size() == 1) {
			if (r1List.size()>1) {
				boolean result = processReadLists(r1List, r2List, tags1, tags2);
				return result;
			} else {
				boolean result = processReadLists(r2List, r1List, tags2, tags1);
				return result;
			}
		}

		// Both lists have multiple reads, sort and compare pairwise
		if (r1List.size() != r2List.size()) {
			throw new IllegalStateException("Both lists have multiple reads, but different sizes.  This should never happen.");
			// return false; // Different sizes cannot match pairwise
		}

		r1List.sort(PAIRED_READ_ORDER_COMPARATOR);
		r2List.sort(PAIRED_READ_ORDER_COMPARATOR);

		for (int i = 0; i < r1List.size(); i++) {
			SAMRecord r1 = r1List.get(i);
			SAMRecord r2 = r2List.get(i);

			// If any pair does not match, return false immediately
			if (!compareTags(r1, r2, tags1, tags2, true)) {
				return false;
			}
		}

		return true; // All pairs matched
	}

	/**
	 * Processes a multi-read list and a single-read list to find concordant reads based on specific tags.
	 *
	 * This method filters the multi-read list based on the presence of all required tags (`tagsMultiread`),
	 * then attempts to find concordant reads between the filtered multi-read list and the single-read list.
	 * If no concordant reads are found, it selects the first record from the filtered multi-read list.
	 *
	 * @param multiReadList   The list of reads containing multiple records.
	 * @param singleReadList  The list of reads containing a single record.
	 * @param tagsMultiread   A list of tags used for filtering and comparison in the multi-read list.
	 * @param tagsSingleRead  A list of tags used for comparison between records in both lists.
	 * @return                A boolean value indicating whether concordant reads were found and processed.
	 */
	public boolean processReadLists(List<SAMRecord> multiReadList, List<SAMRecord> singleReadList,
									List<String> tagsMultiread, List<String> tagsSingleRead) {

		// Filter the multi-read list based on tagsMultiread
		// This gets rid of reads that don't have the tags of interest.
		List<SAMRecord> filteredMultiReadList = multiReadList.stream()
				.filter(record -> tagsMultiread.stream().allMatch(tag -> record.getStringAttribute(tag) != null))
				.collect(Collectors.toList());

		// If the filtered list is empty, add the first record from the original multi-read list
		// if no reads have the tags of interest, then we're going to eventually record that.
		if (filteredMultiReadList.isEmpty()) {
			filteredMultiReadList.add(multiReadList.getFirst());
		}

		// Explicitly handle the case where both lists have exactly one read
		// The data has been simplified back to the base case of a single read in each list that can be compared.
		if (filteredMultiReadList.size() == 1 && singleReadList.size() == 1) {
			return compareTags(filteredMultiReadList.getFirst(), singleReadList.getFirst(), tagsMultiread, tagsSingleRead, true);
		}

		// There's still multiple reads in the multi-read list, so we need to find the concordant read.
		// Attempt to find concordant reads between the filtered multi-read list and the single read
		// if there were multiple reads that could be concordant, try to find the concordant read.
		SAMRecord singleReadRecord = singleReadList.getFirst();
		List<Integer> matchingIndices = IntStream.range(0, filteredMultiReadList.size())
				.filter(index -> compareTags(filteredMultiReadList.get(index), singleReadRecord, tagsMultiread, tagsSingleRead, false))
				.boxed()
				.toList();

		SAMRecord selectedMultiReadRecord;
		if (matchingIndices.isEmpty()) {
			// If no concordant reads, select the first record from the filtered list
			selectedMultiReadRecord = filteredMultiReadList.getFirst();
		} else {
			// If a concordant reads exist, select the first matching record
			selectedMultiReadRecord = filteredMultiReadList.get(matchingIndices.getFirst());
		}

		// Save the tag values for the selected pair of records and return the result of compareTags
		return compareTags(selectedMultiReadRecord, singleReadRecord, tagsMultiread, tagsSingleRead, true);
	}


	private boolean compareTags(final SAMRecord r1, final SAMRecord r2, final List<String> tags1, final List<String> tags2, boolean saveTagValues) {
		List<String> tagValuesList = new ArrayList<>();
		boolean tagValuesConcordant = true;

		for (int i = 0; i < tags1.size(); i++) {
			String v1 = getTagValue(r1, tags1.get(i));
			String v2 = getTagValue(r2, tags2.get(i));

			tagValuesList.add(v1);
			tagValuesList.add(v2);

			// Case 1: Both values are null (concordant)
			if (v1 == null && v2 == null) {
				continue;
			}

			// Case 2: One value is null and the other is not (discordant)
			if (v1 == null || v2 == null) {
				tagValuesConcordant = false;
				continue;
			}

			// Case 3: Both values are non-null; check equality
			if (!v1.equals(v2)) {
				tagValuesConcordant = false;
				continue;
			}
		}

		if (saveTagValues)
			saveTagValues(tagValuesList, tagValuesConcordant);

		return tagValuesConcordant;
	}


	private String getTagValue (final SAMRecord r, final String tag) {
		String result = r.getStringAttribute(tag);
		if (result==null)
			return null;
		if (result.equals(this.TAG_MISSING_VALUE))
			result=null;
		return result;
	}

	private void saveTagValues (List<String> tagValuesList, boolean tagValuesConcordant) {
		if (DISCORDANT_READS_ONLY && tagValuesConcordant)
			return;
		// short circuit if the output is null.
		if (TAG_VALUES_OUTPUT==null)
			return;

		switch (mode) {
			case NO_COMPRESSION -> {
				String tagValues = String.join(TAG_DELIMITER, tagValuesList);
				tagValuesString.increment(tagValues);
			}
			case FIXED_LENGTH_COMPRESSION -> {
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

	/**
	 * Validate that the tag values have consistent lengths so they can be decoded later.
	 * Missing data is encoded with a length of 0 and is always valid.
	 * If a tag length is initialized as non-zero length it will always need to maintain that length.
	 * If the tag length is missing, it will be updated when a non-zero length is observed, but will then be fixed to that length.
	 * @param tagValuesList The list of tag values to validate.
	 * @return True if the tag values have consistent lengths, false otherwise.
	 */
	private boolean validateTagLengths(List<String> tagValuesList) {
		// Compute the lengths of the current tag values, treating nulls as length 0
		int[] lengths = tagValuesList.stream()
				.mapToInt(tag -> tag == null ? 0 : tag.length())
				.toArray();

		// Initialize tagLengths if not already set
		if (this.tagLengths == null) {
			this.tagLengths = lengths;
			return true; // Always valid on first initialization
		}

		// Validate and update tagLengths dynamically
		boolean valid = true;
		for (int i = 0; i < lengths.length; i++) {
			if (lengths[i] == 0) {
				// Current length is 0; no conflict, continue
				continue;
			}

			if (this.tagLengths[i] == 0) {
				// Update stored tagLengths if the current length is non-zero
				this.tagLengths[i] = lengths[i];
			} else if (this.tagLengths[i] != lengths[i]) {
				// Conflict: current length does not match stored value
				valid = false;
			}
		}

		return valid;
	}


	private PeekableIterator<List<SAMRecord>> getReadIterator (SamHeaderAndIterator headerAndIter, final Integer readQuality, List<String> requiredTags) {
		// sort the data by query name.
		// Iterator<SAMRecord> iter = getQueryNameSortedData(headerAndIter);
		Iterator<SAMRecord> iter = headerAndIter.iterator;

		// Filter out reads below a map quality threshold.  This removes non-primary reads.
		iter = new MapQualityFilteredIterator(iter, readQuality, true).iterator();

		// If any tag value is set to the missing value (default "-") remove the read.
		/*
		if (TAG_MISSING_VALUE!=null) {
			for (String tag: requiredTags)
				iter = new TagValueFilteringIterator<String>(iter, tag, List.of(TAG_MISSING_VALUE), false);
		}
		*/

		if (this.REMOVE_TAG_SUFFIX!=null) {
			for (String tag: requiredTags)
				iter = new BAMTagCleanupIterator.Builder(iter)
						.tag(tag)
						.suffixToRemove(REMOVE_TAG_SUFFIX)
						.build();
		}

		//  Don't filter out reads without tags.  Missing tags are considered discordant.
		// 	iter = new MissingTagFilteringIterator(iter, requiredTags.toArray(new String[0]));

		// if there are suffixes on the read names /1 and /2, remove them.
		Pattern pattern = Pattern.compile("/[12]$");
		iter = new ReadNameCleanupIterator.Builder(iter)
				.patternToRemove(pattern)
				.build();

		// Queryname sort data if it is not already sorted by queryname.  This will spill to disk.
		// if the data was already in queryname order, it can be streamed and filtered directly.
		if (!headerAndIter.header.getSortOrder().equals(SAMFileHeader.SortOrder.queryname)) {
			log.info("Input SAM/BAM not in queryname order, sorting...");
			final ProgressLogger progressLogger = new ProgressLogger(log, 1000000, "Sorting reads in query name order");
			iter = SamRecordSortingIteratorFactory.create(headerAndIter.header, iter, READ_NAME_COMPARATOR, progressLogger);
		}

		final GroupingIterator<SAMRecord> groupingIterator = new GroupingIterator<>(iter, READ_NAME_COMPARATOR);
		PeekableIterator<List<SAMRecord>> peekableIterator = new PeekableIterator<>(groupingIterator);
		return peekableIterator;
	}

	/**
	 * A comparator that sorts reads by query name.
	 * @return A metrics file.
	 */
	static final Comparator<SAMRecord> READ_NAME_COMPARATOR =  new Comparator<SAMRecord>() {
		private final SAMRecordQueryNameComparator comp = new SAMRecordQueryNameComparator();
		@Override
		public int compare(final SAMRecord s1, final SAMRecord s2) {
			return comp.fileOrderCompare(s1, s2);
		}
	};

	static final Comparator<SAMRecord> QUERYNAME_COMPARATOR =  new Comparator<SAMRecord>() {
		private final SAMRecordQueryNameComparator comp = new SAMRecordQueryNameComparator();
		@Override
		public int compare(final SAMRecord s1, final SAMRecord s2) {
			return comp.compare(s1, s2);
		}
	};

	/**
	 * A comparator sorts by first of pair/second of pair.
	 */
	static final Comparator<SAMRecord> PAIRED_READ_ORDER_COMPARATOR =  new Comparator<SAMRecord>() {
		private final SAMRecordQueryNameComparator comp = new SAMRecordQueryNameComparator();
		@Override
		public int compare(final SAMRecord s1, final SAMRecord s2) {
			boolean r1Paired = s1.getReadPairedFlag();
			boolean r2Paired = s2.getReadPairedFlag();
			if (r1Paired || r2Paired) {
				if (!r1Paired) {
					return 1;
				}

				if (!r2Paired) {
					return -1;
				}

				if (s1.getFirstOfPairFlag() && s2.getSecondOfPairFlag()) {
					return -1;
				}

				if (s1.getSecondOfPairFlag() && s2.getFirstOfPairFlag()) {
					return 1;
				}
			}
			return 0;
		}
	};



	public enum CompressionMode {
		NO_COMPRESSION(false, false),
		FIXED_LENGTH_COMPRESSION(false, true),
		VARYING_LENGTH_COMPRESSION(true, false);

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

	/**
	 * A custom comparator that sorts concatenated keys by their components.
	 * Use this to sort concatonated tag values by their components consistently to the byte array representation.
	 */
	public static class CustomKeyComparator implements Comparator<String> {
		public static final CustomKeyComparator INSTANCE = new CustomKeyComparator();

		private CustomKeyComparator() {}

		@Override
		public int compare(String key1, String key2) {
			// Split and preprocess the keys
			String[] components1 = splitAndReplaceNulls(key1);
			String[] components2 = splitAndReplaceNulls(key2);

			// Compare component by component
			int len = Math.min(components1.length, components2.length);
			for (int i = 0; i < len; i++) {
				String comp1 = components1[i];
				String comp2 = components2[i];

				// Handle "NA" and "null" values (sort them before any alphanumerics)
				if (comp1.equals("NA") && !comp2.equals("NA")) {
					return -1;
				}
				if (!comp1.equals("NA") && comp2.equals("NA")) {
					return 1;
				}

				// Compare non-"NA" components lexicographically
				int cmp = comp1.compareTo(comp2);
				if (cmp != 0) {
					return cmp;
				}
			}

			// If all compared components are equal, the shorter key is smaller
			return Integer.compare(components1.length, components2.length);
		}

		/**
		 * Splits the concatenated key by ":" and replaces any null or "null" substrings with "NA".
		 *
		 * @param key The concatenated key.
		 * @return An array of substrings with "null" values replaced by "NA".
		 */
		private String[] splitAndReplaceNulls(String key) {
			if (key == null) {
				return new String[]{"NA"}; // Handle completely null key
			}
			return Arrays.stream(key.split(":"))
					.map(component -> component == null || component.equals("null") ? "NA" : component)
					.toArray(String[]::new);
		}
	}



	protected String[] customCommandLineValidation() {

		final ArrayList<String> list = new ArrayList<>(1);


		if (this.TAGS_1.size()!=this.TAGS_2.size())
			list.add("TAGS_1 and TAGS_2 must be the same length.");

		if (BAM_OUTPUT_1 != null)
			IOUtil.assertFileIsWritable(this.BAM_OUTPUT_1);

		if (BAM_OUTPUT_2 != null)
			IOUtil.assertFileIsWritable(this.BAM_OUTPUT_2);

		if (READ_COUNT_OUTPUT != null)
			IOUtil.assertFileIsWritable(this.READ_COUNT_OUTPUT);

		if (TAG_VALUES_OUTPUT != null)
			IOUtil.assertFileIsWritable(this.TAG_VALUES_OUTPUT);

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
