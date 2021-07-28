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

package org.broadinstitute.dropseqrna.sbarro;

import java.io.File;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.Collectors;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.Sbarro;
import org.broadinstitute.dropseqrna.metrics.UmiSharingMetrics;
import org.broadinstitute.dropseqrna.metrics.umisharing.ParentEditDistanceMatcher;
import org.broadinstitute.dropseqrna.metrics.umisharing.ParentEditDistanceMatcher.TagValues;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;
import org.broadinstitute.dropseqrna.utils.MultiComparator;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.ObjectSink;
import org.broadinstitute.dropseqrna.utils.PeekableGroupingIterator;
import org.broadinstitute.dropseqrna.utils.ProgressLoggingIterator;
import org.broadinstitute.dropseqrna.utils.SamHeaderUtil;
import org.broadinstitute.dropseqrna.utils.SamWriterSink;
import org.broadinstitute.dropseqrna.utils.StringInterner;
import org.broadinstitute.dropseqrna.utils.StringTagComparator;
import org.broadinstitute.dropseqrna.utils.editdistance.CollapseTagWithContext;
import org.broadinstitute.dropseqrna.utils.editdistance.FindSimilarEntitiesResult;
import org.broadinstitute.dropseqrna.utils.editdistance.MapBarcodesByEditDistance;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.DefaultTaggingIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.MapQualityPredicate;
import org.broadinstitute.dropseqrna.utils.readiterators.RequiredTagPredicate;
import org.broadinstitute.dropseqrna.utils.readiterators.SamRecordSortingIteratorFactory;

import com.google.common.collect.Sets;

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
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(summary = "Collapse rabies viruses where half the barcode matches within some edit distance, within a cell.",
oneLineSummary = "CollTapse rabies virus barcodes.)",
programGroup = Sbarro.class)
public class BipartiteRabiesVirusCollapse extends CommandLineProgram {

	private static final Log log = Log.getInstance(BipartiteRabiesVirusCollapse.class);
	
	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze.  Must be coordinate sorted. ", optional=false)
	public File INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="A BAM with rabies barcode tags.")
	public File OUTPUT;

	@Argument(doc="A rather verbose report containing the barcode collapse information for each cell.", optional=true)
	public File REPORT;

	@Argument(doc="The cell barcode tag.")
	public String CELL_BARCODE_TAG="XC";

	@Argument (doc="The rabies full barcode tag.  This is the concatonation of the stop codon barcode and polyA barcode.")
	public String RABIES_BARCODE_TAG="rm";

	@Argument(doc="The output tag for the newly collapsed tag values")
	public String OUT_TAG;

	@Argument (doc="Use UMI sharing by itself to determine who should be collapsed.  Ignore similarities (or lack thereof) of rabies virus barcode sequences.")
	public Boolean UMI_SHARING_ONLY_MODE=false;
	
	@Argument(doc="When UMI_SHARING_ONLY_MODE=true, use this umi sharing threshold to determine if barcodes should be collapsed")
	public Double UMI_SHARING_THRESHOLD=0.8;
	
	@Argument(doc="When UMI_SHARING_ONLY_MODE=true, UMIs are considered shared between two rabies viruses if they are within this edit distance.  0 indicates a perfect match.")
	public Integer UMI_SHARING_EDIT_DISTANCE=0;
	
	@Argument (doc="Where to split the rabies barcode for the edit distance check.")
	public Integer RABIES_BARCODE_SPLIT_POSITION=10;
	
	@Argument(doc="The molecular barcode tag.")
	public String MOLECULAR_BARCODE_TAG="XM";

	@Argument(doc="Read quality filter.  Filters all reads lower than this mapping quality.  Defaults to 10.  Set to 0 to not filter reads by map quality.")
	public Integer READ_MQ=10;

	@Argument(doc="The edit distance to collapse rabies viruses, calculating the minimum ED of the halfs of the rabies virus [STOP_CODON_BARCODE_TAG/POLY_A_BARCODE_TAG].")
	public Integer EDIT_DISTANCE=1;

	@Argument (doc="The minumum number of UMIs rabies virus to be eligible for collapse.  Must be >= 1.")
	public int MIN_COUNT=1;

	@Argument(doc = "Should indels be considered in edit distance calculations?  Doing this correctly is far slower than a simple edit distance test, but gives a more complete result.")
	public boolean FIND_INDELS=false;

	@Argument(doc="Number of threads to use.  Defaults to 1.")
	public int NUM_THREADS=1;

	private ForkJoinPool forkJoinPool=null;
	private MapBarcodesByEditDistance med;
	private static DecimalFormat df2 = new DecimalFormat("#.##");

	@Override
	protected int doWork() {
		// Log.setGlobalLogLevel(LogLevel.WARNING);
		if (this.NUM_THREADS>1) this.forkJoinPool = new ForkJoinPool(this.NUM_THREADS);
		IOUtil.assertFileIsReadable(this.INPUT);
		IOUtil.assertFileIsWritable(this.OUTPUT);

		PrintStream outMetrics = null;
		if (this.REPORT!=null) {
			outMetrics = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(this.REPORT));
			if (this.UMI_SHARING_ONLY_MODE) writeMetricsHeaderUMISharing(outMetrics);
			else writeMetricsHeader(outMetrics);
		}
		this.med = new MapBarcodesByEditDistance(false, this.NUM_THREADS, 0);
		IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        SamReader reader = SamReaderFactory.makeDefault().open(INPUT);
        
        SAMFileWriter writer = getWriter (reader);
        final SamWriterSink uninformativeRecSink = new SamWriterSink(writer);
        // final SamReader reader, final String collapseTag, final int mapQuality, String outTag, ObjectSink<SAMRecord> uninformativeReadsSink
        PeekableGroupingIterator<SAMRecord> groupingIter = orderReadsByTagsPeekable(reader, this.RABIES_BARCODE_TAG, this.READ_MQ, this.OUT_TAG, uninformativeRecSink);
                
        log.info("Collapsing tag and writing results");
        // set up filters.
        /*
        int maxNumInformativeReadsInMemory=1000;
        while (groupingIter.hasNext()) {
        	
        	List<SAMRecord> informativeRecs = new ArrayList<>();
        	// you have to grab the next element, in case it's the first of the group but not the first group!
        	informativeRecs.add(groupingIter.next()); 
        	
        	while (groupingIter.hasNextInGroup())         		
        		informativeRecs.add(groupingIter.next()); 
        	
        	// you have all the informative reads.
        	// do some additional logging if the number of reads is bigger than what you've seen before.        	
        	if (informativeRecs.size()>maxNumInformativeReadsInMemory) {
        		maxNumInformativeReadsInMemory=informativeRecs.size();
        		log.info("Max informative reads in memory [" + maxNumInformativeReadsInMemory +"]");        	
        	}
        	
        	// get context.
        	List<SAMRecord> result;
        	if (!this.UMI_SHARING_ONLY_MODE)
        		result = processRecordList(informativeRecs, this.RABIES_BARCODE_TAG, this.OUT_TAG, this.MOLECULAR_BARCODE_TAG, this.EDIT_DISTANCE, this.FIND_INDELS, this.MIN_COUNT, outMetrics);
        	else
        		result = processRecordListUMISharing (informativeRecs, this.RABIES_BARCODE_TAG, this.OUT_TAG, this.MOLECULAR_BARCODE_TAG, this.EDIT_DISTANCE, this.UMI_SHARING_EDIT_DISTANCE, this.MIN_COUNT, outMetrics);
        	// write out the informative reads, then all the reads that were not altered.
        	result.stream().forEach(writer::addAlignment); 
        }
        */
        
        
        while (groupingIter.hasNext()) {
        	SortingCollection<SAMRecord> iter = buildPerContextSortingIterator (groupingIter, reader.getFileHeader());	        	        	        	
        	if (!this.UMI_SHARING_ONLY_MODE)
        		processRecordList(iter, this.RABIES_BARCODE_TAG, this.OUT_TAG, this.MOLECULAR_BARCODE_TAG, this.EDIT_DISTANCE, this.FIND_INDELS, this.MIN_COUNT, outMetrics, writer);
        	else
        		processRecordListUMISharing (iter, this.RABIES_BARCODE_TAG, this.OUT_TAG, this.MOLECULAR_BARCODE_TAG, this.EDIT_DISTANCE, this.UMI_SHARING_EDIT_DISTANCE, this.MIN_COUNT, outMetrics, writer);        	         	
        }
        
        
        log.info("Re-sorting output BAM in genomic order.");
        CloserUtil.close(groupingIter);
        CloserUtil.close(reader);
        writer.close();
        if (outMetrics!=null) CloserUtil.close(outMetrics);
        log.info("DONE");
		return 0;

	}
	
	private SortingCollection<SAMRecord> buildPerContextSortingIterator (PeekableGroupingIterator<SAMRecord> groupingIter, SAMFileHeader header) {
		SortingCollection<SAMRecord> sortingCollection = SortingCollection.newInstance(SAMRecord.class, new BAMRecordCodec(header), CollapseTagWithContext.NO_OP_COMPARATOR, this.MAX_RECORDS_IN_RAM);

    	// you have to grab the next element, in case it's the first of the group but not the first group!
    	sortingCollection.add(groupingIter.next()); 
    	
    	// spool the reads for a whole group into the sorting collection to operate on - the code uses a multi-pass approach so we can't just iterate over the grouping iterator.
    	while (groupingIter.hasNextInGroup())         		
    		sortingCollection.add(groupingIter.next());
    	
    	// wrap up the sorting collection for adding records.
    	sortingCollection.doneAdding();
    	sortingCollection.setDestructiveIteration(false);
		return (sortingCollection);
	}
	
	/**
	 * Process a batch of records and use UMI Sharing alone as the way to decide which rabies genes should be merged.
	 * @param informativeRecs
	 * @param rabiesVirusBarcodeTag
	 * @param outTag
	 * @param molecularBCTag
	 * @param editDistance
	 * @param findIndels
	 * @param minUMICount
	 * @param outMetrics
	 * @return
	 */
	private void processRecordListUMISharing (final Iterable<SAMRecord> i, final String rabiesVirusBarcodeTag, final String outTag, final String molecularBCTag, final int editDistance, final int umiEditDistance, final int minUMICount, final PrintStream out, final SAMFileWriter writer) {
		// make a new peekable iterator once for each time you go through the data and have to start over.
		PeekableIterator<SAMRecord> iter = new PeekableIterator<>(i.iterator());
    	if (!iter.hasNext()) return;
    	    	
    	String cellBarcode = iter.peek().getStringAttribute(this.CELL_BARCODE_TAG);
    	// log.info(cellBarcode);
    	
		ParentEditDistanceMatcher matcher = new ParentEditDistanceMatcher(Collections.singletonList(molecularBCTag), Collections.singletonList(umiEditDistance), false, this.NUM_THREADS);
		double sharingThreshold = 0.8d;
		MapBarcodesByEditDistance mbed = new MapBarcodesByEditDistance(false, this.NUM_THREADS, 10000);	
		Map<String, Set<TagValues>> tagValuesMap = getTagValues(iter, rabiesVirusBarcodeTag, molecularBCTag, matcher);

		// get the count of the tag values.
		ObjectCounter<String> rabiesVirusCounter = new ObjectCounter<>();
		for (String key: tagValuesMap.keySet()) {
			rabiesVirusCounter.incrementByCount(key, tagValuesMap.get(key).size());
		}
		rabiesVirusCounter.filterByMinCount(minUMICount);
		
		FindSimilarEntitiesResult<String, UmiSharingMetrics> result =mbed.collapseBarcodesByUmiSharing(rabiesVirusCounter, matcher, sharingThreshold, tagValuesMap, minUMICount);
		
		// write out metrics
		result.getCollapseMetric().stream().forEach(x -> writeMetricsBody(out, cellBarcode, x));
		
		// remap BAM.
		Map<String, String> map = reverseMap(result.getEntityMap());
		iter = new PeekableIterator<>(i.iterator());
		
		//apply the collapse to the BAM records.
		while (iter.hasNext()) {
			SAMRecord rec = iter.next();
			String oldLabel = rec.getStringAttribute(rabiesVirusBarcodeTag);
			String newLabel = map.get(oldLabel);
			if (newLabel==null) newLabel=oldLabel;
			rec.setAttribute(outTag, newLabel);
			writer.addAlignment(rec);
		}				
	}
	
	private Map<String,String> reverseMap (Map<String, List<String>> map) {
		Map<String,String> result = new HashMap<>();		
		for (String key: map.keySet()) {
			List<String> values = map.get(key);
			if (values!=null && values.size()>0)
				for (String v: values)
					result.put(v, key);				
		}
		return result;
	}
	
	private Map<String, Set<TagValues>> getTagValues (final Iterator<SAMRecord> informativeRecs, final String rabiesVirusBarcodeTag, final String molecularBCTag, ParentEditDistanceMatcher matcher) {
		StringInterner interner = new StringInterner();
		
		Map<String, Set<TagValues>> umis = new HashMap<String, Set<TagValues>>();						
		while (informativeRecs.hasNext()) {
			SAMRecord r = informativeRecs.next();
			String barcode = r.getStringAttribute(rabiesVirusBarcodeTag);
			TagValues v = matcher.getValues(r);
			Set<TagValues> s = umis.get(barcode);
			if (s==null) {
				s = new HashSet<TagValues>();
				barcode=interner.intern(barcode);
				umis.put(barcode, s);
			}			
			s.add(v);
		}
		return umis;
	}
	
	/**
	 * This is the true workhorse of the program.
	 * Try to group rabies virus barcodes based on the edit distance of either of the two portions that make up the rabies virus barcode.
	 * @param informativeRecs
	 * @return
	 */
	private void processRecordList (final Iterable<SAMRecord> informativeRecs, final String rabiesVirusBarcodeTag, final String outTag, final String molecularBCTag, final int editDistance, final boolean findIndels, final int minUMICount, final PrintStream outMetrics, final SAMFileWriter writer) {

		// make a new peekable iterator once for each time you go through the data and have to start over.
		PeekableIterator<SAMRecord> iter = new PeekableIterator<>(informativeRecs.iterator());
		if (!iter.hasNext()) return;
				
    	String cellBarcode = iter.peek().getStringAttribute(this.CELL_BARCODE_TAG);
		Map<String, ObjectCounter<String>> umis= getUMIsForTag (iter, rabiesVirusBarcodeTag, molecularBCTag);
		ObjectCounter<String> rabiesVirusCounter = getCollapedUMICounts(umis, editDistance, minUMICount);
		Map<String, BipartiteRabiesVirusCollapseResultCollection> collapseResult = collapseRabiesBarcodes(rabiesVirusCounter, this.RABIES_BARCODE_SPLIT_POSITION, findIndels, editDistance, this.med);
		// add the umi counts and sharing.
		collapseResult =addUMISharing (collapseResult, umis);

		// write reporting for this cell.
		writeMetricsBody(outMetrics, cellBarcode, collapseResult);

		// get the best mapping.
		Map<String, String> bestMapping = new HashMap<>();
		for (BipartiteRabiesVirusCollapseResultCollection c: collapseResult.values()) {
			BipartiteRabiesVirusCollapseResult best = c.getBestResult();
			if (best!=null)
				bestMapping.put(best.getChildBarcode(), best.getParentBarcode());

		}
		iter = new PeekableIterator<>(informativeRecs.iterator());
		//apply the collapse to the BAM records.
		while (iter.hasNext()) {
			SAMRecord rec = iter.next();
			String oldLabel = rec.getStringAttribute(rabiesVirusBarcodeTag);
			String newLabel = bestMapping.get(oldLabel);
			if (newLabel==null) newLabel=oldLabel;
			rec.setAttribute(outTag, newLabel);
			writer.addAlignment(rec);
		}		
	}
			

	/**
	 *
	 * @param collapseMap
	 * @param umis
	 * @return
	 */
	public Map<String, BipartiteRabiesVirusCollapseResultCollection> addUMISharing (final Map<String, BipartiteRabiesVirusCollapseResultCollection> collapseMap, final Map<String, ObjectCounter<String>> umis) {
		for (String key: collapseMap.keySet()) {
			// need to check for ties and break them (?)
			BipartiteRabiesVirusCollapseResultCollection c  = collapseMap.get(key);
			// double [] umiSharing = c.getResults().stream().mapToDouble((x ->  calculateFractionUMIsShared (x.getChildBarcode(), x.getParentBarcode(), umis))).toArray();
			c.getResults().stream().forEach(x -> x.setUmiSharing(calculateFractionUMIsShared (x.getChildBarcode(), x.getParentBarcode(), umis)));
			c.getResults().stream().forEach(x -> x.setUmisChild(umis.get(x.getChildBarcode()).getSize()));
			c.getResults().stream().forEach(x -> x.setUmisParent(umis.get(x.getParentBarcode()).getSize()));
			BipartiteRabiesVirusCollapseResult best = c.getBestResult();
			// if (best!=null) log.info(best.toString());
		}
		return collapseMap;
	}

	public double calculateFractionUMIsShared (final String barcode1, final String barcode2, final Map<String, ObjectCounter<String>> umis) {
		Set<String> umiSet1 = new HashSet<>(umis.get(barcode1).getKeys());
		Set<String> umiSet2 = new HashSet<>(umis.get(barcode2).getKeys());
		int both= Sets.intersection(umiSet1, umiSet2).size();
		int denom = Math.min(umiSet1.size(), umiSet2.size());
		double r = (double) both / (double) denom;
		return r;
	}

	/**
	 * Map between rabies barcodes that are related on the left or right hand side.
	 *
	 * @param rabiesBarcodes A list of rabies barcodes to find merges in.
	 * @param rabiesSplitPosition A position at which to split the barcodes.
	 * @param findIndels consider indels (false=hamming only)
	 * @param editDistance the edit distance to collapse barcodes at
	 * @param med
	 * @return A map containing a key which is the non-merged barcode, the value is the barcode the key should be merged into.
	 */
	public Map<String, BipartiteRabiesVirusCollapseResultCollection> collapseRabiesBarcodes (final ObjectCounter<String> rabiesBarcodes, final int rabiesSplitPosition, final boolean findIndels, final int editDistance, final MapBarcodesByEditDistance med) {

		// map original sequence to left half.
		Map<String,String> leftHalf = new HashMap<>();
		rabiesBarcodes.getKeysOrderedByCount(true).stream().forEach(x-> leftHalf.put(x, x.substring(0, rabiesSplitPosition)));

		// map original sequence to right half
		Map<String,String> rightHalf = new HashMap<>();
		rabiesBarcodes.getKeysOrderedByCount(true).stream().forEach(x-> rightHalf.put(x, x.substring(rabiesSplitPosition, x.length())));

		List<String> barcodeList = rabiesBarcodes.getKeysOrderedByCount(true);

		// this needs to be a list of results.
		// then we'll need to adjudicate this result.
		// Multimap<String, String> result = ArrayListMultimap.create();

		Map<String, BipartiteRabiesVirusCollapseResultCollection> result = new HashMap<>();

		while (barcodeList.isEmpty()==false) {
			// get the nearby barcodes.
			String b = barcodeList.get(0);
			String left = leftHalf.get(b);
			String right = rightHalf.get(b);

			barcodeList.remove(b);
			// get cached barcode halfs.
			List<String> lefts = barcodeList.stream().map(x -> leftHalf.get(x)).collect(Collectors.toList());
			List<String> rights = barcodeList.stream().map(x -> rightHalf.get(x)).collect(Collectors.toList());

			int [] edLeft = med.getEditDistanceDistributioneMultithreaded(left, lefts, findIndels);
			int [] edRight = med.getEditDistanceDistributioneMultithreaded(right, rights, findIndels);

			// iterate and see if left OR right is close.
			List<String> closeBCList = new ArrayList<>();
			for (int i=0; i<edLeft.length; i++)
				if (edLeft[i]<=editDistance || edRight[i]<=editDistance) {
					String neighbor = barcodeList.get(i);
					BipartiteRabiesVirusCollapseResultCollection c = result.get(neighbor);
					if (c==null) {
						c = new BipartiteRabiesVirusCollapseResultCollection(neighbor);
						result.put(neighbor, c);
					}
					c.add(b, edLeft[i], edRight[i]);
				}
			barcodeList.removeAll(closeBCList);
		}
		return result;
	}

	/**
	 * This gets counts of UMIs for a rabies tag, and does UMI edit distance collapse as part of the counting.
	 * @param informativeRecs
	 * @param tag
	 * @param molBCTag
	 * @param countTagsEditDistance
	 * @param findIndels
	 * @return
	 */
	private Map<String, ObjectCounter<String>> getUMIsForTag (final Iterator<SAMRecord> informativeRecs, final String tag, final String molBCTag) {
		// For each barcode, extract the unique set of countTag values.
		StringInterner interner = new StringInterner();

		// to implement edit distance collapse of tag values, this needs to be an object counter.
		Map<String, ObjectCounter<String>> countTagValues = new HashMap<>();

		while (informativeRecs.hasNext()) {
			SAMRecord r = informativeRecs.next();
			String barcode = r.getStringAttribute(tag);
			ObjectCounter<String> valuesSet = countTagValues.get(barcode);
			// if the set doesn't exist initialize and add...
			if (valuesSet==null) {
				valuesSet=new ObjectCounter<>();
				countTagValues.put(barcode, valuesSet);
			}
			String v = r.getStringAttribute(molBCTag);
			String val = interner.intern(v);
			valuesSet.increment(val);
		}
		return countTagValues;
	}

	private ObjectCounter<String> getCollapedUMICounts (final Map<String, ObjectCounter<String>> countTagValues, final Integer countTagsEditDistance, final Integer minUMICount) {
		// collapse the tag values if needed for each count tag.
		if (countTagsEditDistance>0)
			for (String key: countTagValues.keySet()) {
				ObjectCounter<String> value = countTagValues.get(key);
				value = med.collapseAndMergeBarcodes(value, false, countTagsEditDistance);
				countTagValues.put(key, value);
			}

		// now count the values.
		ObjectCounter<String> barcodeCounts = new ObjectCounter<>();
		for (String barcode: countTagValues.keySet()) {
			// perform collapse here on each object counter.
			int count = countTagValues.get(barcode).getKeys().size();
			barcodeCounts.incrementByCount(barcode, count);
		}
		barcodeCounts.filterByMinCount(minUMICount);
		return barcodeCounts;
	}

	private SAMFileWriter getWriter (final SamReader reader) {
		SAMFileHeader header = reader.getFileHeader();
		SamHeaderUtil.addPgRecord(header, this);

		header.addComment("Edit distance collapsed using rabies virus barcode halves at position [" + this.RABIES_BARCODE_SPLIT_POSITION + "] to new tag " + this.OUT_TAG+ " with edit distance "+ this.EDIT_DISTANCE + "using indels=" + this.FIND_INDELS);
        SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, this.OUTPUT);
        return writer;
	}


	private PeekableGroupingIterator<SAMRecord> orderReadsByTagsPeekable (final SamReader reader, final String collapseTag, final int mapQuality, String outTag, ObjectSink<SAMRecord> uninformativeReadsSink) {
		// SORT on the context tags.
		List <String> contextTags = Arrays.asList(CELL_BARCODE_TAG);
		List<String> expectedTags = Arrays.asList(this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, this.RABIES_BARCODE_TAG);
		
		StringTagComparator [] comparators = contextTags.stream().map(x -> new StringTagComparator(x)).toArray(StringTagComparator[]::new);
		final MultiComparator<SAMRecord> multiComparator = new MultiComparator<>(comparators);

		// set up filters.
        MapQualityPredicate mapQualityPredicate = CollapseTagWithContext.getMapQualityPredicate(mapQuality);
        RequiredTagPredicate requiredTagPredicate = CollapseTagWithContext.getRequiredTagPredicate(collapseTag, expectedTags);
        
        // log progress on read iteration
        ProgressLogger progressLogger = new ProgressLogger(this.log);
        ProgressLoggingIterator progressLoggingIter = new ProgressLoggingIterator(reader.iterator(), progressLogger);
        
        // apply a default result tag to all reads - this is useful for reads that are not in the analysis and automatically sunk to the writer.
        DefaultTaggingIterator iter = new DefaultTaggingIterator(progressLoggingIter.iterator(), collapseTag, outTag);
        
        // reads that don't pass the filter are sunk, reads that pass are sorted and grouped.
		InformativeReadFilter filter = new InformativeReadFilter(iter, uninformativeReadsSink, mapQualityPredicate, requiredTagPredicate);		
		
		CloseableIterator<SAMRecord> sortedIter = SamRecordSortingIteratorFactory.create(
                reader.getFileHeader(), filter.iterator(), multiComparator, new ProgressLogger(log));

		PeekableGroupingIterator<SAMRecord> groupedIterator = new PeekableGroupingIterator<>(sortedIter, multiComparator);
		return groupedIterator;
	}

	private void writeMetricsHeader(final PrintStream outMetrics) {
		String [] header = {"CELL_BARCODE", "CHILD_VIRUS_BARCODE", "PARENT_VIRUS_BARCODE", "LEFT_ED", "RIGHT_ED", "UMI_SHARING", "UMI_CHILD", "UMI_PARENTS", "COLLAPSED_FLAG"};
		outMetrics.println(StringUtils.join(header, "\t"));
	}
		
	private void writeMetricsBody(final PrintStream outMetrics, final String cellBarcode, final Map<String, BipartiteRabiesVirusCollapseResultCollection>collapseResult) {
		for (BipartiteRabiesVirusCollapseResultCollection c: collapseResult.values()) {
			BipartiteRabiesVirusCollapseResult best = c.getBestResult();
			for (BipartiteRabiesVirusCollapseResult result: c.getResults()) {
				List<String> body = Arrays.asList(cellBarcode, result.getChildBarcode(), result.getParentBarcode(), Integer.toString(result.getEditDistanceLeft()),Integer.toString(result.getEditDistanceRight()),
				df2.format(result.getUmiSharing()), Integer.toString(result.getUmisChild()), Integer.toString(result.getUmisParent()), Boolean.toString(result==best));
				outMetrics.println(StringUtils.join(body, "\t"));

			}
		}
	}
	
	private void writeMetricsHeaderUMISharing(final PrintStream outMetrics) {
		String [] header = {"CELL_BARCODE", "PARENT_GENE", "CHILD_GENE", "PARENT_SIZE", "CHILD_SIZE", "NUM_SHARED", "FRAC_SHARED"};
		outMetrics.println(StringUtils.join(header, "\t"));
	}
	
	private void writeMetricsBody (final PrintStream outMetrics, String cellBarcode, UmiSharingMetrics metrics) {
		String [] line = {cellBarcode, metrics.PARENT, metrics.CHILD, Integer.toString(metrics.NUM_PARENT), Integer.toString(metrics.NUM_CHILD), Integer.toString(metrics.NUM_SHARED), String.format("%.3f", metrics.FRAC_SHARED)};
		outMetrics.println(StringUtils.join(line, "\t"));
	}
	
	private class InformativeReadFilter extends FilteredIterator<SAMRecord> {
    	private final MapQualityPredicate mapQualityPredicate;
    	private final RequiredTagPredicate requiredTagPredicate;
    	
		protected InformativeReadFilter(Iterator<SAMRecord> underlyingIterator, ObjectSink<SAMRecord> filteredReadSink, MapQualityPredicate mapQualityPredicate, RequiredTagPredicate requiredTagPredicate ) {
			super(underlyingIterator, filteredReadSink);
			this.mapQualityPredicate=mapQualityPredicate;
			this.requiredTagPredicate=requiredTagPredicate;
		}

		@Override
		public boolean filterOut(SAMRecord rec) {			
			// filter out read if either test fails.
			return (! mapQualityPredicate.test(rec) || !requiredTagPredicate.test(rec));
		} 									    	
    }

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new BipartiteRabiesVirusCollapse().instanceMain(args));
	}


}
