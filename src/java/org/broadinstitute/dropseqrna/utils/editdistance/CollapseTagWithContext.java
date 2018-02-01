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
package org.broadinstitute.dropseqrna.utils.editdistance;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.MultiComparator;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.PeekableGroupingIterator;
import org.broadinstitute.dropseqrna.utils.SamHeaderUtil;
import org.broadinstitute.dropseqrna.utils.StringInterner;
import org.broadinstitute.dropseqrna.utils.StringTagComparator;
import org.broadinstitute.dropseqrna.utils.editdistance.MapBarcodesByEditDistance.AdaptiveMappingResult;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.MapQualityFilteredIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.MissingTagFilteringIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamRecordSortingIteratorFactory;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.StringUtil;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Argument;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(summary = "Collapse set of barcodes that all share the same BAM tags.  For example, collapse all UMIs that have the same cell, gene, and gene strand tags.  This would be equivilent to collapsing the UMIs in DGE.",
oneLineSummary = "Collapse barcodes in the context of one or more tags.)",
programGroup = DropSeq.class)

public class CollapseTagWithContext extends CommandLineProgram {

	private static final Log log = Log.getInstance(CollapseTagWithContext.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.  Must be coordinate sorted. ", optional=false)
	public File INPUT;

	@Argument(doc="Collapse tags that are within <EDIT_DISTANCE>, and have the same CONTEXT_TAGS.  For example, if your context tags were cell and gene, you could collapse UMI tags.", optional=false)
	public String COLLAPSE_TAG;

	@Argument(doc="Group reads by these read tags.  Collapse the COLLAPSE_TAG values that have the same CONTEXT_TAGS values.  Reads with unset CONTEXT_TAGS that will be grouped together and loaded into memory together.  "
			+ "This can cause a large amount of memory usage if you pick a lot of tags that are all mostly not set.", minElements = 1)
	public List<String> CONTEXT_TAGS;

	@Argument (doc="By default, groups of reads are gathered by their CONTEXT_TAGS and ordered by the number of total reads.  Contexts with larger numbers of reads are potential 'parents' of smaller context objects. "
			+ "If this option is used, the count of a context to determine it's ordering is the unique count of values of the TAG(S) added here.  "
			+ "For example, if you wanted to collapse by UMI counts instead of read counts, you could put the UMI tag here.")
	public List<String> COUNT_TAGS;

	@Argument(doc="The output tag for the newly collapsed tag values")
	public String OUT_TAG;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output BAM file with the new collapsed tag.", optional=false)
	public File OUTPUT;

	@Argument(doc="The edit distance to collapse tags.  If adaptive edit distance is used, this is the default edit distance used if no adaptive edit distance is discovered.")
	public Integer EDIT_DISTANCE=1;

	@Argument(doc = "Should indels be considered in edit distance calculations?  Doing this correctly is far slower than a simple edit distance test, but gives a more complete result.")
	public boolean FIND_INDELS=false;

	@Argument(doc="Read quality filter.  Filters all reads lower than this mapping quality.  Defaults to 10.  Set to 0 to not filter reads by map quality.")
	public Integer READ_MQ=10;

	@Argument (doc="The minumum number of reads (unless using the COUNT_TAGS option) for a context to be eligible for collapse.", optional=true)
	public Integer MIN_COUNT=null;

	@Argument(doc="Number of threads to use.  Defaults to 1.")
	public int NUM_THREADS=1;

	@Argument (doc="Instead of using the default fixed edit distance, use an adaptive edit distance.  "
			+ "For each mergable entity, this tries to determine if there are 2 clusters of data by edit distance, and only merge the close-by neighbors.  EXPERIMETNAL!!!")
	public boolean ADAPTIVE_EDIT_DISTANCE=false;

	@Argument (doc="If adaptive edit distance is used, this is the maximum edit distance allowed.", optional=true)
	public Integer ADAPTIVE_ED_MAX=-1;

	@Argument (doc="If adaptive edit distance is used, this is the minimum edit distance allowed.", optional=true)
	public Integer ADAPTIVE_ED_MIN=-1;

	@Argument (doc="If provided, writes out some metrics about each barcode that is merged.", optional=true)
	public File ADAPTIVE_ED_METRICS_FILE;

	// make this once and reuse it.
	private MapBarcodesByEditDistance med;

	@Override
	protected int doWork() {
		if (this.ADAPTIVE_EDIT_DISTANCE & this.ADAPTIVE_ED_MAX==null) {
			log.error("If adaptive edit distance is in use, must set a maximum adaptive edit distance!");
			return 1;
		}
		if (this.ADAPTIVE_EDIT_DISTANCE & this.ADAPTIVE_ED_MIN==null) {
			log.error("If adaptive edit distance is in use, must set a minimum adaptive edit distance!");
			return 1;
		}
		PrintStream outMetrics = null;
		if (this.ADAPTIVE_ED_METRICS_FILE!=null) {
			outMetrics = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(this.ADAPTIVE_ED_METRICS_FILE));
			writeMetricsHeader(outMetrics);
		}


		med = new MapBarcodesByEditDistance(false, this.NUM_THREADS, 0);
		IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        SamReader reader = SamReaderFactory.makeDefault().open(INPUT);
        SAMFileWriter writer = getWriter (reader);

        PeekableGroupingIterator<SAMRecord> groupingIter = orderReadsByTagsPeekable(reader, this.COLLAPSE_TAG, this.CONTEXT_TAGS, this.READ_MQ);
        ProgressLogger pl = new ProgressLogger(log);

        log.info("Collapsing tag and writing results");

        // set up filters.
        FilteredIterator<SAMRecord> mapFilter = getMapFilteringIterator(this.READ_MQ);
        FilteredIterator<SAMRecord> tagFilter = getMissingTagIterator(this.COLLAPSE_TAG, this.CONTEXT_TAGS);

        int maxNumInformativeReadsInMemory=1000;

        while (groupingIter.hasNext()) {

        	// prime the first read of the group.
        	SAMRecord r = groupingIter.next();
        	List<SAMRecord> informativeRecs = new ArrayList<>();
        	informativeRecs = getInformativeRead(r, informativeRecs, this.COLLAPSE_TAG, this.OUT_TAG, mapFilter, tagFilter, writer, pl);

        	// get subsequent reads.  Informative reads stay in memory, uninformative reads are written to disk and discarded.
        	while (groupingIter.hasNextInGroup()) {
        		r = groupingIter.next();
        		informativeRecs = getInformativeRead(r, informativeRecs, this.COLLAPSE_TAG, this.OUT_TAG, mapFilter, tagFilter, writer, pl);
        	}

        	// you have all the informative reads.
        	// do some additional logging if the number of reads is bigger than what you've seen before.
        	boolean verbose = false;
        	if (informativeRecs.size()>maxNumInformativeReadsInMemory) {
        		maxNumInformativeReadsInMemory=informativeRecs.size();
        		log.info("Max informative reads in memory [" + maxNumInformativeReadsInMemory +"]");
        		verbose=true;
        	}

        	// remove the informative reads, as they will get retagged.
        	List<SAMRecord> result = processRecordList(informativeRecs, this.COLLAPSE_TAG, this.COUNT_TAGS, this.OUT_TAG, this.FIND_INDELS, this.EDIT_DISTANCE, this.ADAPTIVE_ED_MIN, this.ADAPTIVE_ED_MAX, this.MIN_COUNT, verbose, outMetrics);

        	// write out the informative reads, then all the reads that were not altered.
        	result.stream().forEach(writer::addAlignment);
        }
        log.info("Re-sorting output BAM in genomic order.");
        CloserUtil.close(groupingIter);
        CloserUtil.close(reader);
        writer.close();
        if (outMetrics!=null) CloserUtil.close(outMetrics);
        log.info("DONE");
		return 0;
	}

	/**
	 * Tests to see if a read is informative.
	 * If it is, add it to the list of informative reads.
	 * If not, write the record to the writer and don't alter the informative reads list.
	 * @param r
	 * @param mapFilter
	 * @param tagFilter
	 * @param writer
	 * @return
	 */
	private List<SAMRecord> getInformativeRead (final SAMRecord r, final List<SAMRecord> informativeReads, final String collapseTag, final String outTag,
			final FilteredIterator<SAMRecord> mapFilter, final FilteredIterator<SAMRecord> tagFilter, final SAMFileWriter writer, final ProgressLogger pl) {
		pl.record(r);
		boolean informative = !mapFilter.filterOut(r) && !tagFilter.filterOut(r);
		if (!informative) {
			// if the uninformative read has a value for the collapse tag, set the output tag to the same value.
			String v = r.getStringAttribute(collapseTag);
			if (v!=null)
				r.setAttribute(outTag, v);
			writer.addAlignment(r);
		}
		else
			informativeReads.add(r);
		return informativeReads;
	}

	private List<SAMRecord> processRecordList (final Collection<SAMRecord> informativeRecs, final String collapseTag, final List<String> countTags, final String outTag, final boolean findIndels, final int editDistance, final Integer minEditDistance, final Integer maxEditDistance, final Integer minNumObservations, final boolean verbose, final PrintStream outMetrics) {
		if (informativeRecs.size()==0) return Collections.EMPTY_LIST;

		// get context.
		String context = getContextString(informativeRecs.iterator().next(), this.CONTEXT_TAGS);

		// get barcode counts.
		ObjectCounter<String> barcodeCounts = getBarcodeCounts (informativeRecs, collapseTag, countTags);
		Map<String, String> collapseMap = collapseBarcodes(barcodeCounts, findIndels, editDistance, minEditDistance, maxEditDistance, verbose, outMetrics, context, minNumObservations);

		// now that you have a map from children to the parent, retag all the reads.
		List<SAMRecord> result = new ArrayList<>();
		for (SAMRecord r: informativeRecs) {
			String tagValue = r.getStringAttribute(collapseTag);
			// if the tag was not set, then don't set it.
			if (tagValue!=null) {
				// tag was set, set it.
				// is it in the map?  If so, update the tag value.
				if (collapseMap.containsKey(tagValue))
					tagValue = collapseMap.get(tagValue);
				r.setAttribute(outTag, tagValue);
			}
			result.add(r);
		}
		return result;
	}

	private ObjectCounter<String> getBarcodeCounts (final Collection<SAMRecord> informativeRecs, final String collapseTag, final List<String> countTags) {
		// collapse barcodes based on informative reads that have the necessary tags.
		// this counts 1 per read.
		if (countTags==null) {
			List<String> barcodes = informativeRecs.stream().map(x -> x.getStringAttribute(collapseTag)).collect(Collectors.toList());
			ObjectCounter<String> barcodeCounts = new ObjectCounter<>();
			barcodes.stream().forEach(x -> barcodeCounts.increment(x));
			return barcodeCounts;
		}

		// otherwise, for each barcode, extract the unique set of countTag values.
		StringInterner interner = new StringInterner();

		Map<String, Set<String>> countTagValues = new HashMap<>();

		for (SAMRecord r: informativeRecs) {
			String barcode = r.getStringAttribute(collapseTag);
			Set<String> valuesSet = countTagValues.get(barcode);
			// if the set doesn't exist initialize and add...
			if (valuesSet==null) {
				valuesSet=new HashSet<>();
				countTagValues.put(barcode, valuesSet);
			}
			List<String> valsList = new ArrayList<>();
			for (String countTag: countTags) {
				String v = r.getStringAttribute(countTag);
				if (v!=null) valsList.add(v);
			}
			String val = interner.intern(StringUtils.join(valsList, ":"));
			valuesSet.add(val);
		}
		// now count the values.
		ObjectCounter<String> barcodeCounts = new ObjectCounter<>();
		for (String barcode: countTagValues.keySet()) {
			int count = countTagValues.get(barcode).size();
			barcodeCounts.incrementByCount(barcode, count);
		}
		return barcodeCounts;
	}

	private SAMFileWriter getWriter (final SamReader reader) {
		SAMFileHeader header = reader.getFileHeader();
		SamHeaderUtil.addPgRecord(header, this);
		String context = StringUtil.join(" ", this.CONTEXT_TAGS);
		header.addComment("Edit distance collapsed tag " +  this.COLLAPSE_TAG + " to new tag " + this.OUT_TAG+ " with edit distance "+ this.EDIT_DISTANCE + "using indels=" + this.FIND_INDELS + " in the context of tags [" + context + "]");
        SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, this.OUTPUT);
        return writer;
	}

	private Map<String, String> collapseBarcodes(final ObjectCounter<String> barcodeCounts, final boolean findIndels, final int editDistance, final Integer minEditDistance, final Integer maxEditDistance, final boolean verbose, final PrintStream outMetrics, final String context, final Integer minNumObservations) {
		// order the barcodes by the number of reads each barcode has.
		if (verbose) log.info("Collapsing [" + barcodeCounts.getSize() +"] barcodes.");
		if (minNumObservations!=null) barcodeCounts.filterByMinCount(minNumObservations);

		// map of primary barcode to list of child barcodes.
		Map<String, String> result = new HashMap<>();
		Map<String, List<String>> r = null;
		if (ADAPTIVE_EDIT_DISTANCE) {
			AdaptiveMappingResult amr= med.collapseBarcodesAdaptive(barcodeCounts, findIndels, editDistance, minEditDistance, maxEditDistance);
			r = amr.getBarcodeCollapseResult();
			writeMetrics(context, amr, outMetrics);
		}
		else r = med.collapseBarcodes(barcodeCounts, findIndels, editDistance);

		// flip map to each child barcode that belongs to a parent.
		for (String key: r.keySet())
			for (String value: r.get(key))
				result.put(value, key);

		return (result);
	}

	private String getContextString (final SAMRecord r, final List<String> contextTags) {
		List<String> result = new ArrayList<>();
		for (String c: contextTags) {
			String v = r.getStringAttribute(c);
			result.add(v);
		}
		return StringUtils.join(result, ",");
	}

	private void writeMetricsHeader (final PrintStream out) {
		String [] header = {"CONTEXT", "COLLAPSE", "NUM_COLLAPSED", "ADAPTIVE_ED_DISCOVERED", "ADAPTIVE_ED_USED", "NUM_OBS_ORIGINAL", "NUM_OBS_MERGED"};
		out.println(StringUtil.join("\t", header));
	}

	private void writeMetrics (final String context, final AdaptiveMappingResult r, final PrintStream out) {
		if (out==null) return;
		List<EditDistanceMappingMetric> metricList= r.getMetricResult();

		for (EditDistanceMappingMetric edmm: metricList) {
			edmm.getOriginalObservations();
			// Steve reports the number of barcodes including the one that everything is merged into.
			String [] line = {context, edmm.getBarcode(), Integer.toString(edmm.getNumMergedBarcodes()+1), Integer.toString(edmm.getEditDistanceDiscovered()), Integer.toString(edmm.getEditDistanceUsed()),
					Integer.toString(edmm.getOriginalObservations()), Integer.toString(edmm.getTotalObservations())};
			out.println(StringUtil.join("\t", line));
		}
	}

	/**
	 * Given a list of tags, order the reads by the tags, then group by the tags.
	 * This does NOT filter reads.
	 * @param reader
	 * @param tags
	 * @param mapQuality
	 * @return
	 */
	private GroupingIterator<SAMRecord> orderReadsByTags (final SamReader reader, final String collapseTag, final List<String> contextTag, final int mapQuality) {
		// SORT on the context tags.
		StringTagComparator [] comparators = contextTag.stream().map(x -> new StringTagComparator(x)).toArray(StringTagComparator[]::new);
		final MultiComparator<SAMRecord> multiComparator = new MultiComparator<>(comparators);

		CloseableIterator<SAMRecord> sortedIter = SamRecordSortingIteratorFactory.create(
                reader.getFileHeader(), reader.iterator(), multiComparator, new ProgressLogger(log));

		GroupingIterator<SAMRecord> groupedIterator = new GroupingIterator<>(sortedIter, multiComparator);
		return groupedIterator;
	}

	private PeekableGroupingIterator<SAMRecord> orderReadsByTagsPeekable (final SamReader reader, final String collapseTag, final List<String> contextTag, final int mapQuality) {
		// SORT on the context tags.
		StringTagComparator [] comparators = contextTag.stream().map(x -> new StringTagComparator(x)).toArray(StringTagComparator[]::new);
		final MultiComparator<SAMRecord> multiComparator = new MultiComparator<>(comparators);

		CloseableIterator<SAMRecord> sortedIter = SamRecordSortingIteratorFactory.create(
                reader.getFileHeader(), reader.iterator(), multiComparator, new ProgressLogger(log));

		PeekableGroupingIterator<SAMRecord> groupedIterator = new PeekableGroupingIterator<>(sortedIter, multiComparator);
		return groupedIterator;
	}

	private boolean testReadInformative (final SAMRecord r, final FilteredIterator<SAMRecord> mapFilter, final FilteredIterator<SAMRecord> tagFilter) {
		return (!mapFilter.filterOut(r) && !tagFilter.filterOut(r));
	}

	private FilteredIterator<SAMRecord> getMapFilteringIterator (final int mapQuality) {
		FilteredIterator<SAMRecord> mapFilteringIterator = new MapQualityFilteredIterator(Collections.emptyIterator(), mapQuality, false);
		return mapFilteringIterator;
	}

	private FilteredIterator<SAMRecord> getMissingTagIterator (final String collapseTag, final List<String> contextTag) {
		List<String> allTags = new ArrayList<>(contextTag);
		allTags.add(collapseTag);
		String[] tagArray = allTags.stream().toArray(String[]::new);
		FilteredIterator<SAMRecord> missingTagIterator = new MissingTagFilteringIterator(Collections.emptyIterator(), tagArray);
		return missingTagIterator;
	}


	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new CollapseTagWithContext().instanceMain(args));
	}

}
