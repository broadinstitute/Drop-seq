package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.BAMRecordCodec;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMProgramRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;

import java.util.List;

public class CustomBAMIterators {
	
	private static final Log log = Log.getInstance(CustomBAMIterators.class);
	public static final int MAX_RECORDS_IN_RAM = 500000;
	
	public static CloseableIterator<SAMRecord> getReadsInTagOrder (SamReader reader, String primaryTag) {
		// SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
		SAMSequenceDictionary dict= reader.getFileHeader().getSequenceDictionary();
		List<SAMProgramRecord> programs =reader.getFileHeader().getProgramRecords();
		
		final SAMFileHeader writerHeader = new SAMFileHeader();
        writerHeader.setSortOrder(SAMFileHeader.SortOrder.queryname);
        writerHeader.setSequenceDictionary(dict);
        for (SAMProgramRecord spr : programs) {
        	writerHeader.addProgramRecord(spr);
        }
		SortingCollection<SAMRecord> alignmentSorter = SortingCollection.newInstance(SAMRecord.class,
	            new BAMRecordCodec(writerHeader), new BAMTagComparator(primaryTag),
	                MAX_RECORDS_IN_RAM);
		
		log.info("Reading in records for TAG name sorting");
		int counter=0;
		for (SAMRecord r: reader) {
			alignmentSorter.add(r);
			counter++;
			if (counter%1000000==0) log.info(counter+ " records processed");
		}
	
		// CloserUtil.close(reader);
		CloseableIterator<SAMRecord> result = alignmentSorter.iterator();
		log.info("Sorting finished.");
		return (result);
	}
	
	public static CloseableIterator<SAMRecord> getReadsInTagOrder (SamReader reader, List<String> tags) {
		// SamReader reader = SamReaderFactory.makeDefault().open(bamFile);
		
		SAMSequenceDictionary dict= reader.getFileHeader().getSequenceDictionary();
		List<SAMProgramRecord> programs =reader.getFileHeader().getProgramRecords();
		
		final SAMFileHeader writerHeader = new SAMFileHeader();
        writerHeader.setSortOrder(SAMFileHeader.SortOrder.queryname);
        writerHeader.setSequenceDictionary(dict);
        for (SAMProgramRecord spr : programs) {
        	writerHeader.addProgramRecord(spr);
        }
		SortingCollection<SAMRecord> alignmentSorter = SortingCollection.newInstance(SAMRecord.class,
	            new BAMRecordCodec(writerHeader), new BAMTagComparator(tags),
	                MAX_RECORDS_IN_RAM);
		
		log.info("Reading in records for TAG name sorting");
		
		ProgressLogger prog = new ProgressLogger(log);
		
		for (SAMRecord r: reader) {
			alignmentSorter.add(r);
			prog.record(r);
			
		}
	
		CloserUtil.close(reader);
		CloseableIterator<SAMRecord> result = alignmentSorter.iterator();
		log.info("Sorting finished.");
		return (result);
	} 
	
	/**
	 * If the file is sorter in query name order, return an iterator over
	 * the file.  Otherwise, sort records in queryname order and return an iterator over those records.
	 * @param bamFile The BAM/SAM file to read
	 * @return An iterator over the records in the file, in queryname order.
	 */
	public static CloseableIterator<SAMRecord> getQuerynameSortedRecords(SamReader reader) {
		CloseableIterator<SAMRecord> iter = null;
		if (reader.getFileHeader().getSortOrder().equals(SortOrder.queryname)) {
			iter = reader.iterator();
			return iter;
		}
		log.info("Input SAM/BAM not in queryname order, sorting...");
		ProgressLogger p = new ProgressLogger(log, 1000000, "Sorting reads in query name order");
		SAMSequenceDictionary dict= reader.getFileHeader().getSequenceDictionary();
		List<SAMProgramRecord> programs =reader.getFileHeader().getProgramRecords();
		
		final SAMFileHeader writerHeader = new SAMFileHeader();
        writerHeader.setSortOrder(SAMFileHeader.SortOrder.queryname);
        writerHeader.setSequenceDictionary(dict);
        for (SAMProgramRecord spr : programs) {
        	writerHeader.addProgramRecord(spr);
        }
        SortingCollection<SAMRecord> alignmentSorter = SortingCollection.newInstance(SAMRecord.class,
                    new BAMRecordCodec(writerHeader), new SAMRecordQueryNameComparator(),
                        MAX_RECORDS_IN_RAM);
        log.info("Reading in records for query name sorting");
        int counter=0;
        for (SAMRecord r: reader) {
        	p.record(r);
        	alignmentSorter.add(r);
        	counter++;
        	if (counter%1000000==0) log.info(counter+ " records processed");
        }
        CloserUtil.close(reader);
        CloseableIterator<SAMRecord> result = alignmentSorter.iterator();
        log.info("Sorting finished.");
        return result; 
	}
	
}
