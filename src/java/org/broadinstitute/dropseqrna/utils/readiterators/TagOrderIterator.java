package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriterImpl;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.broadinstitute.dropseqrna.utils.bamtagcomparator.BAMTagComparator;
import org.broadinstitute.dropseqrna.utils.bamtagcomparator.ComparatorAggregator;

/**
 * Given an input BAM file, a list of one or more tags, and a list of 0 or more processors
 * @author nemesh
 *
 */
public class TagOrderIterator implements CloseableIterator<SAMRecord> {

	private static final Log log = Log.getInstance(TagOrderIterator.class);	
	private static final ProgressLogger prog = new ProgressLogger(log);
	private final PeekableIterator<SAMRecord> iter;
	private final List<Short> sortingTags;
	private final List<Short> requiredTags;
	
	/**
	 * Set up an iterator that returns reads sorting by tags.
	 * @param inFile the input BAM/SAM file to iterate over
	 * @param sortingTags A list of one of more String value tags to sort the BAM file by
	 * @param requiredTags A list of zero or more String value tags that must have a value set on the read 
	 * @param filters A collection of filters that can filter or add reads to the iterator
	 * @param skipReadsWithoutTags skip all reads that don't have a value set for all of the sorting tags.
	 */
	public TagOrderIterator(File inFile, List<String> sortingTags, List<String> requiredTags, ComparatorAggregator ag, ReadProcessorCollection filters, boolean skipReadsWithoutTags) {
		int maxReads=SAMFileWriterImpl.getDefaultMaxRecordsInRam();
		this.sortingTags=DEIteratorUtils.getShortBAMTags(sortingTags);
		this.requiredTags=DEIteratorUtils.getShortBAMTags(requiredTags);
		this.iter = new PeekableIterator<SAMRecord>(getReadsInTagOrder(inFile, this.sortingTags, this.requiredTags, ag, filters, maxReads, skipReadsWithoutTags));
	}

	/**
	 * Set up an iterator that returns reads sorting by tags.
	 * @param inFile the input BAM/SAM file to iterate over
	 * @param sortingTags A list of one of more String value tags to sort the BAM file by
	 * @param filters A single filter that can filter or add reads to the iterator
	 * @param skipReadsWithoutTags skip all reads that don't have a value set for all of the sorting tags.
	 */
	public TagOrderIterator (File inFile, List<String> sortingTags, List<String> requiredTags, ComparatorAggregator ag, SAMReadProcessorI filter, boolean skipReadsWithoutTags) {
		int maxReads=SAMFileWriterImpl.getDefaultMaxRecordsInRam();
		this.sortingTags=DEIteratorUtils.getShortBAMTags(sortingTags);
		this.requiredTags=DEIteratorUtils.getShortBAMTags(requiredTags);
		ReadProcessorCollection filters = new ReadProcessorCollection();
		filters.addFilter (filter);
		this.iter = new PeekableIterator<SAMRecord>(getReadsInTagOrder(inFile, this.sortingTags, this.requiredTags, ag, filters, maxReads, skipReadsWithoutTags));
	}
	
	/**
	 * Get the BAM tags this iterator is sorted on, in Short format. 
	 * @return
	 */
	public List<Short> getShortTags() {
		return this.sortingTags;
	}
	
	/**
	 * Get the BAM tags this iterator is sorted on, in String format
	 * @return
	 */
	public List<String> getTags () {
		return DEIteratorUtils.getStringBAMTags(this.sortingTags);
	}
	
	/**
	 * This sets up the guts of the object, the CloseableIterator that will be queried by next or hasNext.
	 * @param inFile
	 * @param sortingTags
	 * @param filters
	 * @param maxReadsInRAM
	 * @return
	 */
	private static CloseableIterator<SAMRecord> getReadsInTagOrder (File inFile, List<Short> sortingTags, List<Short> requiredTags, ComparatorAggregator ag, ReadProcessorCollection filters, int maxReadsInRAM, boolean skipReadsWithoutTags) {	
		SamReader reader = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE).open(inFile);
		
		SAMSequenceDictionary dict= reader.getFileHeader().getSequenceDictionary();
		final SAMFileHeader writerHeader = new SAMFileHeader();
        writerHeader.setSequenceDictionary(dict);

        SortingCollection<SAMRecord> alignmentSorter = SortingCollection.newInstance(SAMRecord.class,
	            new BAMRecordCodec(writerHeader), new BAMTagComparator(sortingTags, ag),
	            maxReadsInRAM);
		
		log.info("Reading in records for TAG name sorting");
				
		Collection <SAMRecord> tempList = new ArrayList<SAMRecord>(10);
		int numReadsAdded=0;
		int newNumReads=0;
		for (SAMRecord r: reader) {
			prog.record(r);
			// see if the read has all the sortings tags set and reads without them should be skipped. 
			if (skipReadsWithoutTags  && !testAllAttributesSet(r, requiredTags)) continue;
			tempList  = filters.processRead(r);
			if (tempList.size()>0) numReadsAdded++;
			newNumReads+=tempList.size();
			for (SAMRecord rr: tempList) {
				alignmentSorter.add(rr);
			}
		}
		
		log.info("From BAM file added " + numReadsAdded + " reads out of " +prog.getCount() + " to iterator");
		log.info("New iterator has " + newNumReads + " to iterate on.");
		if (numReadsAdded==0) log.warn("The number of reads added to the iterator was 0");
		CloserUtil.close(reader);
		CloseableIterator<SAMRecord> result = alignmentSorter.iterator();
		log.info("Sorting finished.");
		return (result);
	}
	
	
	/**
	 * For a record and a list of sorting tags, test to see that all BAM tags have a non-null value.
	 * @param r The SAMRecord to test
	 * @param sortingTags a list of BAM tags to test for non-null values.
	 * @return
	 */
	private static boolean testAllAttributesSet (SAMRecord r, List<Short> sortingTags) {
		for (Short s: sortingTags) {
			Object o = r.getAttribute(s);
			if (o==null) return false;
		}
		return true;
	}
	
	
	@Override
	public boolean hasNext() {
		return this.iter.hasNext();
	}

	@Override
	public SAMRecord next() {
		return this.iter.next();		
	}

	@Override
	public void remove() {
		throw new UnsupportedOperationException("Remove not supported");
	}

	@Override
	public void close() {
		this.iter.close();
	}
	
	public SAMRecord peek () {
		return this.iter.peek();
	}

	
	
}
