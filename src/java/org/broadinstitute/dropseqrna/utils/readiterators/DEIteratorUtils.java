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
package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMFileWriterImpl;
import htsjdk.samtools.SAMTagUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;

import java.util.ArrayList;
import java.util.List;

public class DEIteratorUtils {

	private static final Log log = Log.getInstance(DEIteratorUtils.class);
	private ProgressLogger progress = new ProgressLogger(log, 100000);
	
	private Integer MAX_RECORDS_IN_RAM = SAMFileWriterImpl.getDefaultMaxRecordsInRam();
	
	
	
	public static List<Short> getShortBAMTags (List<String> tags) {
		List<Short> result = new ArrayList<Short>(tags.size());
		
		for (String tag : tags) {
			short s = SAMTagUtil.getSingleton().makeBinaryTag(tag);
			result.add(s);
		}
		return (result);
	}
	
	public static List<String> getStringBAMTags (List<Short> tags) {
		List<String> result = new ArrayList<String>(tags.size());
		
		for (Short tag : tags) {
			String s = SAMTagUtil.getSingleton().makeStringTag(tag);
			result.add(s);
		}
		return (result);
	}
	
	
	
	
	/**
	 * Utility method to get a collection of records from a peekable iterator.
	 * @param iter
	 * @param geneExonTag
	 * @param cellBarcodeTag
	 * @return
	 */
	/*
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
	
	
	
	
	
	public CloseableIterator<SAMRecord> getReadsInTagOrder (File inFile, List<Short> sortingTags, SAMReadProcessorI filter) {
		ReadProcessorCollection c = new ReadProcessorCollection();		
		c.addFilter (filter);
		return (getReadsInTagOrder(inFile, sortingTags, c));
	}
	
	public CloseableIterator<SAMRecord> getReadsInTagOrder (File inFile, List<Short> sortingTags, ReadProcessorCollection filters) {	
		SamReader reader = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE).open(inFile);
		
		SAMSequenceDictionary dict= reader.getFileHeader().getSequenceDictionary();
		//List<SAMProgramRecord> programs =reader.getFileHeader().getProgramRecords();
		
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
	*/
	
	
	
	/**
	 * When you sort an iterator on multiple fields and some of the fields [1 or more] may be null, this skips past all the null entries and starts with the first non-null entry.
	 * @param iter The iterator to...iterate on
	 * @param emptyAttribute The attribute(s) of the read to scan.
	 * 
	 * @return
	 */
	/*
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
	
	*/

	
		
}
	
