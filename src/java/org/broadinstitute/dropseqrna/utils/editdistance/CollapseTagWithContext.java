package org.broadinstitute.dropseqrna.utils.editdistance;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.MultiComparator;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.StringTagComparator;
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
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(usage = "Collapse set of barcodes that all share the same BAM tags.  For example, collapse all UMIs that have the same cell, gene, and gene strand tags.  This would be equivilent to collapsing the UMIs in DGE.",
usageShort = "Collapse barcodes in the context of one or more tags.)",
programGroup = DropSeq.class)

public class CollapseTagWithContext extends CommandLineProgram {

	private static final Log log = Log.getInstance(CollapseTagWithContext.class);

	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.  Must be coordinate sorted. ", optional=false)
	public File INPUT;

	@Option(doc="Collapse tags that are within <EDIT_DISTANCE>, and have the same CONTEXT_TAGS.  For example, if your context tags were cell and gene, you could collapse UMI tags.", optional=false)
	public String COLLAPSE_TAG;

	@Option(doc="Group reads by these read tags.  Collapse the COLLAPSE_TAG values that have the same CONTEXT_TAGS values.  Reads with unset CONTEXT_TAGS that will be grouped together and loaded into memory together.  "
			+ "This can cause a large amount of memory usage if you pick a lot of tags that are all mostly not set.", minElements = 1)
	public List<String> CONTEXT_TAGS;

	@Option(doc="The output tag for the newly collapsed tag values")
	public String OUT_TAG;

	@Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output BAM file with the new collapsed tag.", optional=false)
	public File OUTPUT;

	@Option(doc="The edit distance to collapse tags")
	public Integer EDIT_DISTANCE=1;

	@Option(doc = "Should indels be considered in edit distance calculations?  Doing this correctly is far slower than a simple edit distance test, but gives a more complete result.")
	public boolean FIND_INDELS=false;

	@Option(doc="Read quality filter.  Filters all reads lower than this mapping quality.  Defaults to 10.  Set to 0 to not filter reads by map quality.")
	public Integer READ_MQ=10;

	@Option(doc="Number of threads to use.  Defaults to 1.")
	public int NUM_THREADS=1;

	private CollapseBarcodeThreaded cbt=null;
	private int threadedBlockSize=20000;

	// make this once and reuse it.
	private MapBarcodesByEditDistance med;

	@Override
	protected int doWork() {
		med = new MapBarcodesByEditDistance(false, this.NUM_THREADS, 0);

		// TODO Auto-generated method stub
		IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        SamReader reader = SamReaderFactory.makeDefault().open(INPUT);
        SAMFileWriter writer = getWriter (reader);

        GroupingIterator<SAMRecord> groupingIter = orderReadsByTags(reader, this.COLLAPSE_TAG, this.CONTEXT_TAGS, this.READ_MQ);
        Iterator<List<SAMRecord>> iter = groupingIter.iterator();
        ProgressLogger pl = new ProgressLogger(log);

        log.info("Collapsing tag and writing results");

        // set up filters.
        FilteredIterator<SAMRecord> mapFilter = getMapFilteringIterator(this.READ_MQ);
        FilteredIterator<SAMRecord> tagFilter = getMissingTagIterator(this.COLLAPSE_TAG, this.CONTEXT_TAGS);

        while (iter.hasNext()) {
        	// get all reads.  This can be semi-big if tags are passed in that aren't set.
        	List<SAMRecord> allRecs = iter.next();
        	// informative reads have all of their tags set, so they can participate in collapse.
        	List<SAMRecord> informativeRecs = allRecs.stream().filter(x -> !mapFilter.filterOut(x)).filter(x -> !tagFilter.filterOut(x)).collect(Collectors.toList());

        	// remove the informative reads, as they will get retagged.
        	allRecs.removeAll(informativeRecs);
        	List<SAMRecord> result = processRecordList(informativeRecs, this.COLLAPSE_TAG, this.OUT_TAG, this.FIND_INDELS, this.EDIT_DISTANCE);

        	// write out the informative reads, then all the reads that were not altered.
        	result.stream().forEach(writer::addAlignment);
        	pl.record(result.stream().toArray(SAMRecord[]::new));
        	allRecs.stream().forEach(writer::addAlignment);
        	pl.record(allRecs.stream().toArray(SAMRecord[]::new));
        }
        log.info("Re-sorting output BAM in genomic order.");
        CloserUtil.close(groupingIter);
        CloserUtil.close(reader);
        writer.close();
        log.info("DONE");
		return 0;
	}

	private List<SAMRecord> processRecordList (final List<SAMRecord> informativeRecs, final String collapseTag, final String outTag, final boolean findIndels, final int editDistance) {

		// collapse barcodes based on informative reads that have the necessary tags.
		List<String> barcodes = informativeRecs.stream().map(x -> x.getStringAttribute(collapseTag)).collect(Collectors.toList());
		Map<String, String> collapseMap = collapseBarcodes(barcodes, findIndels, editDistance);

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

	private SAMFileWriter getWriter (final SamReader reader) {
		SAMFileHeader header = reader.getFileHeader();
		String context = StringUtil.join(" ", this.CONTEXT_TAGS);
		header.addComment("Edit distance collapsed tag " +  this.COLLAPSE_TAG + " to new tag " + this.OUT_TAG+ " with edit distance "+ this.EDIT_DISTANCE + "using indels=" + this.FIND_INDELS + " in the context of tags [" + context + "]");
        SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, this.OUTPUT);
        return writer;
	}

	private Map<String, String> collapseBarcodes(final List<String> barcodes, final boolean findIndels, final int editDistance) {
		// order the barcodes by the number of reads each barcode has.
		ObjectCounter<String> barcodeCounts = new ObjectCounter<>();
		barcodes.stream().forEach(x -> barcodeCounts.increment(x));

		Map<String, String> result = new HashMap<>();


		// map of primary barcode to list of child barcodes.
		Map<String, List<String>> r = med.collapseBarcodes(barcodeCounts, findIndels, editDistance);
		// flip map to each child barcode that belongs to a parent.
		for (String key: r.keySet())
			for (String value: r.get(key))
				result.put(value, key);
		return (result);
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
