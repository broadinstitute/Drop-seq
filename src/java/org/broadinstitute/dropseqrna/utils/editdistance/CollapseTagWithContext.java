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
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.DigitalExpression.DESummary;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;
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

import htsjdk.samtools.BAMRecordCodec;
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
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;
import picard.cmdline.CommandLineProgram;
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

	@Argument (doc="If COUNT_TAGS is set and COUNT_TAGS_EDIT_DISTANCE>0, then collapse the COUNT_TAGS in a CONTEXT by the given edit distance.  For example, if you wanted to collapse "
			+ "by UMIs instead of read counts, and you wanted to further collapse UMIs by edit distance 1, you'd set COUNT_TAGS_EDIT_DISTANCE to 1.  This doesn't do much unless MIN_COUNT is also set "
			+ "as collapse would only be affected if there is a minimum number of counts for each CONTEXT to be in a COLLAPSE.", optional=true)
	public Integer COUNT_TAGS_EDIT_DISTANCE=0;

	@Argument(doc="The output tag for the newly collapsed tag values")
	public String OUT_TAG;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output BAM file with the new collapsed tag.", optional=false)
	public File OUTPUT;

	@Argument(doc="The edit distance to collapse tags.  If adaptive edit distance is used, this is the default edit distance used if no adaptive edit distance is discovered.  If mutational collapse is used, this is the maximum edit distance two barcodes in a network can be apart (but they must have network neighbors at ED=1 for the entire path).")
	public Integer EDIT_DISTANCE=1;

	@Argument(doc = "Should indels be considered in edit distance calculations?  Doing this correctly is far slower than a simple edit distance test, but is a more aggressive method that may be useful in some situations.")
	public boolean FIND_INDELS=false;

	@Argument(doc="Read quality filter.  Filters all reads lower than this mapping quality.  Defaults to 10.  Set to 0 to not filter reads by map quality.")
	public Integer READ_MQ=10;

	@Argument (doc="The minumum number of reads (unless using the COUNT_TAGS option) for a context to be eligible for collapse.  Must be >= 1.")
	public int MIN_COUNT=1;

	@Argument (doc="When collapsing a CONTEXT_TAG, do not emit CONTEXT reads that have fewer than MIN_COUNT counts.  "
			+ "For example, if your context tags were cell and gene and you were collapsing UMI tags and had a MIN_COUNT of 5, then cell/gene pairs with fewer than 5 UMIs "
			+ "would not have their reads emiited in the output BAM.", optional=false)
	public Boolean DROP_SMALL_COUNTS=false;

	@Argument(doc="Number of threads to use.  Defaults to 1.")
	public int NUM_THREADS=1;

	@Argument (doc="Instead of using the default fixed edit distance, use an adaptive edit distance.  "
			+ "For each mergable entity, this tries to determine if there are 2 clusters of data by edit distance, and only merge the close-by neighbors.")
	public boolean ADAPTIVE_EDIT_DISTANCE=false;
	
	@Argument (doc="If adaptive edit distance is used, this is the maximum edit distance allowed.", optional=true)
	public Integer ADAPTIVE_ED_MAX=-1;

	@Argument (doc="If adaptive edit distance is used, this is the minimum edit distance allowed.", optional=true)
	public Integer ADAPTIVE_ED_MIN=-1;

	@Argument (doc="If provided, writes out some metrics about each barcode that is merged by adaptive edit distance collapse.", optional=true)
	public File ADAPTIVE_ED_METRICS_FILE;

	@Argument (doc="If true, add an additional column that contains a comma separated list of edit distances from the current CONTEXT_TAG to all other CONTEXT_TAGS.  This will make files significantly larger!")
	public boolean ADAPTIVE_ED_METRICS_ED_LIST=false;

	@Argument (doc="Instead of using the default fixed edit distance, use a mutational collapse strategy.  "
			+ "For the single largest barcode in the context, find all neighbors within ED=1.  Then find neighbors to those neighbors at ED=1 that are ALSO ED=2 to the original barcode.  Search out to a maximum edit distance of EDIT_DISTANCE.")
	public boolean MUTATIONAL_COLLAPSE=false;
	
	@Argument (doc="If provided, writes out some metrics about each barcode that is merged by mutational edit distance collapse.", optional=true)
	public File MUTATIONAL_COLLAPSE_METRICS_FILE;
	
	@Argument (doc="Use less memory but more time.  Useful if your context groups are huge - very large cells with lots of sequence data, etc.")
	public Boolean LOW_MEMORY_MODE=false;

	// make this once and reuse it.
	private MapBarcodesByEditDistance med;
	private MapBarcodesByEditDistance medUMI;

	int validateCommands () {
		if (this.ADAPTIVE_EDIT_DISTANCE & this.ADAPTIVE_ED_MAX==null) {
			log.error("If adaptive edit distance is in use, must set a maximum adaptive edit distance!");
			return 1;
		}
		if (this.ADAPTIVE_EDIT_DISTANCE & this.ADAPTIVE_ED_MIN==null) {
			log.error("If adaptive edit distance is in use, must set a minimum adaptive edit distance!");
			return 1;
		}
		if (this.MIN_COUNT < 1) {
			log.error(String.format("MIN_COUNT(%d) < 1 does not make sense.", MIN_COUNT));
			return 1;
		}
		if (this.DROP_SMALL_COUNTS && (this.MIN_COUNT < 2)) {
			log.error("If DROP_SMALL_COUNTS is set to true, must set a MIN_COUNT VALUE greater than 1.");
			return 1;
		}
		if (this.COUNT_TAGS_EDIT_DISTANCE>0 && this.COUNT_TAGS==null) {
			log.error ("Edit distance for COUNT_TAGS set, but no COUNT TAGS set.  Can't do edit distance collapse on read counts!");
			return 1;
		}
		if (this.MUTATIONAL_COLLAPSE && this.ADAPTIVE_EDIT_DISTANCE) {
			log.error("Can't specifiy both adaptive edit distance collapse AND mutational collapse.");
			return 1;
		}
		return 0;
	}
	@Override
	protected int doWork() {
		int vc = validateCommands();
		if (vc>0) return vc;

		if (this.COUNT_TAGS_EDIT_DISTANCE>0) this.medUMI = new MapBarcodesByEditDistance(false);

		med = new MapBarcodesByEditDistance(false, this.NUM_THREADS, 0);
		
		PrintStream outMetrics = null;
		if (this.ADAPTIVE_ED_METRICS_FILE!=null) {
			outMetrics = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(this.ADAPTIVE_ED_METRICS_FILE));
			writeAdaptiveEditDistanceMetricsHeader(this.ADAPTIVE_ED_METRICS_ED_LIST, outMetrics);
		}
		
		if (this.MUTATIONAL_COLLAPSE_METRICS_FILE!=null) {
			med = new MapBarcodesByEditDistance(true, this.NUM_THREADS, 1000);
			outMetrics = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(this.MUTATIONAL_COLLAPSE_METRICS_FILE));
			writeMutationalCollapseMetricsHeader(this.ADAPTIVE_ED_METRICS_ED_LIST, outMetrics);
		}

		
		IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        SamReader reader = SamReaderFactory.makeDefault().open(INPUT);
        SAMFileHeader header =  reader.getFileHeader();
        SAMFileWriter writer = getWriter (reader);

        PeekableGroupingIterator<SAMRecord> groupingIter = orderReadsByTagsPeekable(reader, this.COLLAPSE_TAG, this.CONTEXT_TAGS, this.READ_MQ);        

        log.info("Collapsing tag and writing results");

        // set up filters.
        FilteredIterator<SAMRecord> mapFilter = getMapFilteringIterator(this.READ_MQ);
        FilteredIterator<SAMRecord> tagFilter = getMissingTagIterator(this.COLLAPSE_TAG, this.CONTEXT_TAGS);

        if (!LOW_MEMORY_MODE) 
        	fasterIteration(groupingIter, mapFilter, tagFilter, writer, outMetrics);
        else
        	lowMemoryIteration(groupingIter, mapFilter, tagFilter, writer, outMetrics, header);
        
        
        log.info("Re-sorting output BAM in genomic order.");
        CloserUtil.close(groupingIter);
        CloserUtil.close(reader);
        writer.close();
        if (outMetrics!=null) CloserUtil.close(outMetrics);
        log.info("DONE");
		return 0;
	}
	
	/**
	 * With this method, we keep all the records in memory 
	 * @param groupingIter
	 * @param mapFilter
	 * @param tagFilter
	 * @param writer
	 * @param outMetrics
	 */
	private void fasterIteration (PeekableGroupingIterator<SAMRecord> groupingIter, FilteredIterator<SAMRecord> mapFilter, FilteredIterator<SAMRecord> tagFilter,
			SAMFileWriter writer, PrintStream outMetrics) {
		log.info("Running fast single iteration mode");
		int maxNumInformativeReadsInMemory=1000;
		ProgressLogger pl = new ProgressLogger(log);
        while (groupingIter.hasNext()) {
        	
        	// prime the first read of the group.
        	SAMRecord r = groupingIter.next();
        	List<SAMRecord> informativeRecs = new ArrayList<>();
        	r = getInformativeRead(r, this.COLLAPSE_TAG, this.OUT_TAG, mapFilter, tagFilter, writer, pl);
        	if (r!=null) informativeRecs.add(r);
        	// get subsequent reads.  Informative reads stay in memory, uninformative reads are written to disk as part of the final BAM output and discarded.
        	while (groupingIter.hasNextInGroup()) {
        		r = groupingIter.next();
        		r = getInformativeRead(r, this.COLLAPSE_TAG, this.OUT_TAG, mapFilter, tagFilter, writer, pl);
        		if (r!=null) informativeRecs.add(r);
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
        	processRecordList(informativeRecs, writer, this.COLLAPSE_TAG, this.COUNT_TAGS, this.OUT_TAG,
					this.FIND_INDELS, this.EDIT_DISTANCE, this.ADAPTIVE_ED_MIN, this.ADAPTIVE_ED_MAX, this.MIN_COUNT,
                    this.DROP_SMALL_COUNTS, this.COUNT_TAGS_EDIT_DISTANCE, verbose, outMetrics, this.ADAPTIVE_ED_METRICS_ED_LIST);
        	
        }		
	}
	
	private void lowMemoryIteration (PeekableGroupingIterator<SAMRecord> groupingIter, FilteredIterator<SAMRecord> mapFilter, FilteredIterator<SAMRecord> tagFilter,
			SAMFileWriter writer, PrintStream outMetrics, SAMFileHeader header) {
		log.info("Running (slower) memory efficient mode");
		ProgressLogger pl = new ProgressLogger(log);
		boolean verbose = false;
        while (groupingIter.hasNext()) {   
        	// for this group, get a SortingCollection.
        	SortingCollection<SAMRecord> sortingCollection = SortingCollection.newInstance(SAMRecord.class, new BAMRecordCodec(header), NO_OP_COMPARATOR, this.MAX_RECORDS_IN_RAM);
        	
        	// prime the first read of the group.
        	SAMRecord r = groupingIter.next();
        	r = getInformativeRead(r, this.COLLAPSE_TAG, this.OUT_TAG, mapFilter, tagFilter, writer, pl);        	
        	if (r!=null) sortingCollection.add(r);
        	
        	while (groupingIter.hasNextInGroup()) {
        		r = groupingIter.next();
        		r = getInformativeRead(r, this.COLLAPSE_TAG, this.OUT_TAG, mapFilter, tagFilter, writer, pl);
        		if (r!=null) sortingCollection.add(r);
        	}
        	// wrap up the sorting collection for adding records.
        	sortingCollection.doneAdding();
        	sortingCollection.setDestructiveIteration(false);
        	
        	PeekableIterator<SAMRecord> iter = new PeekableIterator<>(sortingCollection.iterator());
        	if (!iter.hasNext()) return;

    		// get context.
    		String context = getContextString(iter.peek(), this.CONTEXT_TAGS);

    		// get barcode counts.
    		ObjectCounter<String> barcodeCounts = getBarcodeCounts (iter, this.COLLAPSE_TAG, this.COUNT_TAGS, this.COUNT_TAGS_EDIT_DISTANCE);
    		if (this.MIN_COUNT > 1 & !this.MUTATIONAL_COLLAPSE) barcodeCounts.filterByMinCount(this.MIN_COUNT);
    		
    		Map<String, String> collapseMap = collapseBarcodes(barcodeCounts, this.FIND_INDELS, this.EDIT_DISTANCE, this.ADAPTIVE_ED_MIN, this.ADAPTIVE_ED_MAX, this.MIN_COUNT, verbose, outMetrics, context, this.ADAPTIVE_ED_METRICS_ED_LIST);
    		iter = new PeekableIterator<>(sortingCollection.iterator());
    		retagBarcodedReads(iter, barcodeCounts, collapseMap, this.DROP_SMALL_COUNTS, writer, this.COLLAPSE_TAG, this.OUT_TAG);        	
        	
        }
		
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
	public static SAMRecord getInformativeRead (final SAMRecord r, final String collapseTag, final String outTag,
			final FilteredIterator<SAMRecord> mapFilter, final FilteredIterator<SAMRecord> tagFilter, final SAMFileWriter writer, final ProgressLogger pl) {
		pl.record(r);
		boolean informative = !mapFilter.filterOut(r) && !tagFilter.filterOut(r);
		if (!informative) {
			// if the uninformative read has a value for the collapse tag, set the output tag to the same value.
			String v = r.getStringAttribute(collapseTag);
			if (v!=null)
				r.setAttribute(outTag, v);
			writer.addAlignment(r);
			return null;
		}					
		return r;
	}

	// Break this up into 2 methods - one to collapse results, one to retag reads.
	// Make the input a peekable iterator instead of a collection of SAMRecords.
	private void processRecordList (final Collection<SAMRecord> informativeRecs, final SAMFileWriter writer, final String collapseTag,
                                               final List<String> countTags, final String outTag, final boolean findIndels,
                                               final int editDistance, final Integer minEditDistance, final Integer maxEditDistance,
                                               final int minNumObservations, final boolean dropSmallCounts,
                                               final Integer countTagsEditDistance, final boolean verbose, final PrintStream outMetrics, final boolean writeEditDistanceDistribution) {
		if (informativeRecs.size()==0) return;

		// get context.
		String context = getContextString(informativeRecs.iterator().next(), this.CONTEXT_TAGS);

		// get barcode counts.
		ObjectCounter<String> barcodeCounts = getBarcodeCounts (new PeekableIterator<SAMRecord>(informativeRecs.iterator()), collapseTag, countTags, countTagsEditDistance);
		if (minNumObservations > 1 & !this.MUTATIONAL_COLLAPSE) barcodeCounts.filterByMinCount(minNumObservations);
		
		Map<String, String> collapseMap = collapseBarcodes(barcodeCounts, findIndels, editDistance, minEditDistance, maxEditDistance, minNumObservations, verbose, outMetrics, context, writeEditDistanceDistribution);
		retagBarcodedReads(informativeRecs.iterator(), barcodeCounts, collapseMap, dropSmallCounts, writer, collapseTag, outTag);
	}
	
	private void retagBarcodedReads (Iterator<SAMRecord> informativeRecs, ObjectCounter<String> barcodeCounts, Map<String, String> collapseMap, boolean dropSmallCounts, SAMFileWriter writer,
			String collapseTag, String outTag) {
		
		Set<String> expectedBarcodes = null;
		// already validated that if dropSmallCounts is true, then the minNumObservations > 1.
		if (dropSmallCounts)
			// use all the remaining barcodes that have counts.
			expectedBarcodes = new HashSet<>(barcodeCounts.getKeys());
		while (informativeRecs.hasNext()) {
			SAMRecord r = informativeRecs.next();
			String tagValue = r.getStringAttribute(collapseTag);
			// if the tag was not set, then don't set it.
			if (tagValue!=null) {
				// tag was set,
				// if the tagValue is not in the expected barocde list and the barcode list is populated, then don't add this read and short circuit to next read.
				if (expectedBarcodes!=null && !expectedBarcodes.contains(tagValue))
					continue;
				// tag was set, set it.
				// is it in the map?  If so, update the tag value.
				if (collapseMap.containsKey(tagValue))
					tagValue = collapseMap.get(tagValue);
				r.setAttribute(outTag, tagValue);
			}
			writer.addAlignment(r);
		}		
	}


	/**
	 * Operate on a collection of records, find the counts of each context.
	 * This may include fancy work like edit distance collapse of the count tags.
	 * @param informativeRecs
	 * @param collapseTag
	 * @param countTags
	 * @param countTagsEditDistance
	 * @return
	 */
	private ObjectCounter<String> getBarcodeCounts (final PeekableIterator<SAMRecord> informativeRecs, final String collapseTag, final List<String> countTags, final Integer countTagsEditDistance) {
		// collapse barcodes based on informative reads that have the necessary tags.
		// this counts 1 per read.
		if (countTags.isEmpty()) {
			List<String> barcodes = informativeRecs.stream().map(x -> x.getStringAttribute(collapseTag)).collect(Collectors.toList());
			ObjectCounter<String> barcodeCounts = new ObjectCounter<>();
			barcodes.stream().forEach(x -> barcodeCounts.increment(x));
			return barcodeCounts;
		}

		// otherwise, for each barcode, extract the unique set of countTag values.
		StringInterner interner = new StringInterner();

		// to implement edit distance collapse of tag values, this needs to be an object counter.
		Map<String, ObjectCounter<String>> countTagValues = new HashMap<>();
		
		while (informativeRecs.hasNext()) {
			SAMRecord r=informativeRecs.next();
			String barcode = r.getStringAttribute(collapseTag);
			ObjectCounter<String> valuesSet = countTagValues.get(barcode);
			// if the set doesn't exist initialize and add...
			if (valuesSet==null) {
				valuesSet=new ObjectCounter<>();
				countTagValues.put(barcode, valuesSet);
			}
			// if there are multiple count tags, need to distinguish between them.  IE: if your count was of distinct UMI + some strand tag, then you'd need a distinct list of those 2 tags aggregated together, and the count
			// is the number of unique values.
			List<String> valsList = new ArrayList<>();
			for (String countTag: countTags) {
				String v = r.getStringAttribute(countTag);
				if (v!=null) valsList.add(v);
			}
			String val = interner.intern(StringUtils.join(valsList, ":"));

			valuesSet.increment(val);
		}

		// collapse the tag values if needed for each count tag.
		if (countTagsEditDistance>0)
			for (String key: countTagValues.keySet()) {
				ObjectCounter<String> value = countTagValues.get(key);
				value = medUMI.collapseAndMergeBarcodes(value, false, countTagsEditDistance);
				countTagValues.put(key, value);
			}

		// now count the values.
		ObjectCounter<String> barcodeCounts = new ObjectCounter<>();
		for (String barcode: countTagValues.keySet()) {
			// perform collapse here on each object counter.
			int count = countTagValues.get(barcode).getKeys().size();
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

	private Map<String, String> collapseBarcodes(final ObjectCounter<String> barcodeCounts, final boolean findIndels, final Integer editDistance, final Integer minEditDistance, final Integer maxEditDistance, final Integer minSizeToCollapse, final boolean verbose, final PrintStream outMetrics, final String context, final boolean writeEditDistanceDistribution) {
		// order the barcodes by the number of reads each barcode has.
		if (verbose) log.info("Collapsing [" + barcodeCounts.getSize() +"] barcodes.");

		// map of primary barcode to list of child barcodes.
		Map<String, String> result = new HashMap<>();
		Map<String, List<String>> r = null;
		if (this.ADAPTIVE_EDIT_DISTANCE && !this.MUTATIONAL_COLLAPSE) {
			AdaptiveMappingResult amr= med.collapseBarcodesAdaptive(barcodeCounts, findIndels, editDistance, minEditDistance, maxEditDistance);
			r = amr.getBarcodeCollapseResult();
			writeMetrics(writeEditDistanceDistribution, context, amr, outMetrics);
		} else if (this.MUTATIONAL_COLLAPSE && !this.ADAPTIVE_EDIT_DISTANCE) {
			r=med.collapseBarcodesByMutationalCollapse(barcodeCounts, findIndels, editDistance, minSizeToCollapse);
			ObjectCounter<String> aggregatedCounts=aggregateCounts(barcodeCounts, r);
			writeMutationalReport(barcodeCounts, aggregatedCounts, r, outMetrics);
		} 
		else r = med.collapseBarcodes(barcodeCounts, findIndels, editDistance);
		

		// flip map to each child barcode that belongs to a parent.
		for (String key: r.keySet())
			for (String value: r.get(key))
				result.put(value, key);

		return (result);
	}
	
	/**
	 * Generate a new object counter reflects the counts of barcodes after collapsing those barcodes via a given mapping.
	 * @param data The counts of barcodes (UMIs/reads/etc).
	 * @param mapping The map of a barcode (key) to some other set of barcodes (values) that collapse into that barcode.
	 * @return The counts of barcodes after collapse - some barcodes with have higher counts, some barcodes that were collapsed will be removed from the original data.  
	 * The sum of total counts will be the same as the input data.
	 */
	public static ObjectCounter<String> aggregateCounts (ObjectCounter<String> data, Map<String, List<String>> mapping) {
		ObjectCounter<String> result = new ObjectCounter<>();
		// build out all the new counts for barcodes that have merged results.
		Set<String> mergedBarcodes=new HashSet<String>();
		
		for (String key: mapping.keySet()) {
			List<String> values = mapping.get(key);
			int totalCount=data.getCountForKey(key);			
			mergedBarcodes.addAll(values);
			totalCount += values.stream().mapToInt(x-> data.getCountForKey(x)).sum();						
			result.incrementByCount(key, totalCount);
			mergedBarcodes.add(key);			
		}
		
		// build out the singletons that were not merged by finding the non-merged data keys.
		for (String k: data.getKeys()) 
			if (!mergedBarcodes.contains(k)) {
				result.incrementByCount(k, data.getCountForKey(k));
			}
		
		return result;		
	}

	private String getContextString (final SAMRecord r, final List<String> contextTags) {
		List<String> result = new ArrayList<>();
		for (String c: contextTags) {
			String v = r.getStringAttribute(c);
			result.add(v);
		}
		return StringUtils.join(result, ",");
	}

	private void writeAdaptiveEditDistanceMetricsHeader (final boolean writeEditDistanceDistribution, final PrintStream out) {
		List<String> header = new ArrayList<String>(Arrays.asList("CONTEXT", "COLLAPSE", "NUM_COLLAPSED", "ADAPTIVE_ED_DISCOVERED", "ADAPTIVE_ED_USED", "NUM_OBS_ORIGINAL", "NUM_OBS_MERGED"));
		if (writeEditDistanceDistribution)
			header.add("ED_DISTRIBUTION");
		out.println(StringUtil.join("\t", header));
	}

	private void writeMutationalCollapseMetricsHeader (final boolean writeEditDistanceDistribution, final PrintStream out) {
		String [] header = {"sequence",  "counts", "parent",  "edist",  "fam_seqs", "fam_counts"};
		out.println(StringUtils.join(header, "\t"));
	}

	private void writeMetrics (final boolean writeEditDistanceDistribution, final String context, final AdaptiveMappingResult r, final PrintStream out) {
		if (out==null) return;
		List<EditDistanceMappingMetric> metricList= r.getMetricResult();

		for (EditDistanceMappingMetric edmm: metricList) {
			edmm.getOriginalObservations();
			// Steve reports the number of barcodes including the one that everything is merged into.
			List<String> line = new ArrayList<>(Arrays.asList(context, edmm.getBarcode(), Integer.toString(edmm.getNumMergedBarcodes()+1), Integer.toString(edmm.getEditDistanceDiscovered()), Integer.toString(edmm.getEditDistanceUsed()),
					Integer.toString(edmm.getOriginalObservations()), Integer.toString(edmm.getTotalObservations())));

			if (writeEditDistanceDistribution) {
				int [] edList = edmm.getEdList();
				if (edList.length>0) {
					Integer[] x = Arrays.stream( edList ).boxed().toArray( Integer[]::new );
					String edFormatted = StringUtil.join(",", x);
					line.add(edFormatted);
				} else
					line.add("NA");

			}
			out.println(StringUtil.join("\t", line));

		}
	}
	
	private void writeMutationalReport (ObjectCounter<String> data, ObjectCounter<String> aggregateCounts, Map<String, List<String>> mapping, PrintStream out) {						
		for (String parentSeq: mapping.keySet()) {		
			List<String> sequences = mapping.get(parentSeq);														
			int famSeqs=sequences.size()+1;
			String [] line = {parentSeq, Integer.toString(data.getCountForKey(parentSeq)), parentSeq, "0", Integer.toString(famSeqs) ,Integer.toString(aggregateCounts.getCountForKey(parentSeq))};
			out.println(StringUtils.join(line, "\t"));
															
			for (String v: sequences) {
				int ed = HammingDistance.getHammingDistance(parentSeq, v);
				// for merged results, the family seqs size is always 1 and the fam counts is always 0.
				String [] line2 = {v, Integer.toString(data.getCountForKey(v)), parentSeq, Integer.toString(ed), "1", "0"};
				out.println(StringUtils.join(line2, "\t"));																				
			}			
		}						
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

	public static boolean testReadInformative (final SAMRecord r, final FilteredIterator<SAMRecord> mapFilter, final FilteredIterator<SAMRecord> tagFilter) {
		return (!mapFilter.filterOut(r) && !tagFilter.filterOut(r));
	}

	public static FilteredIterator<SAMRecord> getMapFilteringIterator (final int mapQuality) {
		FilteredIterator<SAMRecord> mapFilteringIterator = new MapQualityFilteredIterator(Collections.emptyIterator(), mapQuality, false);
		return mapFilteringIterator;
	}

	public static FilteredIterator<SAMRecord> getMissingTagIterator (final String collapseTag, final List<String> contextTag) {
		List<String> allTags = new ArrayList<>(contextTag);
		allTags.add(collapseTag);
		String[] tagArray = allTags.stream().toArray(String[]::new);
		FilteredIterator<SAMRecord> missingTagIterator = new MissingTagFilteringIterator(Collections.emptyIterator(), tagArray);
		return missingTagIterator;
	}
	
	static final Comparator<SAMRecord> NO_OP_COMPARATOR =  new Comparator<SAMRecord>() {
        @Override
		public int compare(final SAMRecord e1, final SAMRecord e2) {
            return 0;
        }
    };


	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new CollapseTagWithContext().instanceMain(args));
	}

}
