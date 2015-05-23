package org.broadinstitute.dropseqrna.utils.readIterators;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriterImpl;
import htsjdk.samtools.SAMProgramRecord;
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
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.broadinstitute.dropseqrna.barnyard.Utils;
import org.broadinstitute.dropseqrna.utils.BAMTagComparator;

public class DEIteratorUtils {

	private static final Log log = Log.getInstance(DEIteratorUtils.class);
	private ProgressLogger progress = new ProgressLogger(log, 100000);
	
	private Integer MAX_RECORDS_IN_RAM;
	
	public DEIteratorUtils () {
		MAX_RECORDS_IN_RAM = SAMFileWriterImpl.getDefaultMaxRecordsInRam();
	}
	
	public DEIteratorUtils (int maxRecordsInRam) {
		MAX_RECORDS_IN_RAM = maxRecordsInRam;
	}
	
	
	/**
	 * Utility method to get a collection of records from a peekable iterator.
	 * @param iter
	 * @param geneExonTag
	 * @param cellBarcodeTag
	 * @return
	 */
	
	public Collection<SAMRecord> getRecordCollection (PeekableIterator<SAMRecord> iter, List<String> tags) {
		if (iter.hasNext()==false) return (null);
		List<SAMRecord> result = new ArrayList<SAMRecord>();
		
		SAMRecord r = iter.peek();
		
		List<String> currentValues = getValuesForTags(tags, r);
		while (iter.hasNext()) {
			r=iter.peek();
			
			List<String> nextValues = getValuesForTags(tags, r);
			if (testTagsNotEqual(currentValues, nextValues)) {
				break;
			}
			// this is the same set of records as before, keep going.
			// grab this record for "real" so peek gets the next record that might be in the same gene.
			iter.next();
			this.progress.record(r);
			result.add(r);
			
		}
		return (result);		
	}
	
	private List<String> getValuesForTags(List<String>tags, SAMRecord r) {
		List<String> currentValues = new ArrayList<String>();
		for (String t: tags) {
			currentValues.add(r.getStringAttribute(t));
		}
		return (currentValues);
	}
	
	private boolean testTagsNotEqual (List<String> original, List<String> next) {
		for (int i=0; i<original.size(); i++) {
			String s1 = original.get(i);
			String s2 = next.get(i);
			if (!s1.equals(s2)) return (true);
		}
		return (false);
	}
	
	public CloseableIterator<SAMRecord> getReadsInTagOrder (File inFile, List<String> sortingTags, SAMReadProcessorI filter) {
		ReadProcessorCollection c = new ReadProcessorCollection();		
		c.addFilter (filter);
		return (getReadsInTagOrder(inFile, sortingTags, c));
	}
	
	public CloseableIterator<SAMRecord> getReadsInTagOrder (File inFile, List<String> sortingTags, ReadProcessorCollection filters) {	
		SamReader reader = SamReaderFactory.makeDefault().open(inFile);
		
		SAMSequenceDictionary dict= reader.getFileHeader().getSequenceDictionary();
		List<SAMProgramRecord> programs =reader.getFileHeader().getProgramRecords();
		
		final SAMFileHeader writerHeader = new SAMFileHeader();
        // writerHeader.setSortOrder(SAMFileHeader.SortOrder.queryname);
        writerHeader.setSequenceDictionary(dict);
        for (SAMProgramRecord spr : programs) {
        	writerHeader.addProgramRecord(spr);
        }
		SortingCollection<SAMRecord> alignmentSorter = SortingCollection.newInstance(SAMRecord.class,
	            new BAMRecordCodec(writerHeader), new BAMTagComparator(sortingTags),
	                MAX_RECORDS_IN_RAM);
		
		log.info("Reading in records for TAG name sorting");
		
		ProgressLogger prog = new ProgressLogger(log);
		
		Collection <SAMRecord> tempList = new ArrayList<SAMRecord>(10);
		int numReadsAdded=0;
		for (SAMRecord r: reader) {
			prog.record(r);
			tempList  = filters.processRead(r);
			if (tempList.size()>0) numReadsAdded++;
			for (SAMRecord rr: tempList) {
				alignmentSorter.add(rr);
			}
		}
		log.info("Added " + numReadsAdded + " to iterator out of " +prog.getCount());
		
		if (numReadsAdded==0) log.warn("The number of reads added to the iterator was 0.  Either you have all low map quality data, or there's no geneExon or cellBarcode tags on the reads");
		CloserUtil.close(reader);
		CloseableIterator<SAMRecord> result = alignmentSorter.iterator();
		log.info("Sorting finished.");
		return (result);
	}
	
	
	
	
	/**
	 * When you sort an iterator on multiple fields and some of the fields [1 or more] may be null, this skips past all the null entries and starts with the first non-null entry.
	 * @param iter The iterator to...iterate on
	 * @param emptyAttribute The attribute(s) of the read to scan.
	 * 
	 * @return
	 */
	public PeekableIterator<SAMRecord> primeIterator (PeekableIterator<SAMRecord> iter, String...emptyAttribute) {
		ProgressLogger primeLog = new ProgressLogger(log, 1000000, "Skipped records without tags "+ getFormattedString(emptyAttribute));
		if (iter.hasNext()==false) return (iter);
		
		SAMRecord r = iter.peek();
		// seek to the first gene.
		while (iter.hasNext()) {
			r=iter.peek();
			int numNotNull=0;
			for (String key: emptyAttribute) {
				Object value = r.getAttribute(key);
				if (value!=null) numNotNull++;
			}
			if (numNotNull==emptyAttribute.length){ 
				break;
			}
			r=iter.next();
			primeLog.record(r);
		}
		return (iter);
	}
	
	public PeekableIterator<SAMRecord> primeIterator (PeekableIterator<SAMRecord> iter, List<String> emptyAttributeList) {
		String[] emptyAttribute = emptyAttributeList.toArray(new String[emptyAttributeList.size()]);
		return (primeIterator(iter, emptyAttribute));
	}
	
	private String getFormattedString (String...x) {
		StringBuilder b= new StringBuilder();
		for (int i=0; i<x.length; i++) {
			b.append(x[i]);
			if (i<x.length) {
				b.append(",");
			}
		}
		return (b.toString());
	}
	
	

	
		
}
	