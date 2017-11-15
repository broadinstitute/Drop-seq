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
package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.*;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.util.*;
import org.broadinstitute.dropseqrna.utils.readiterators.SamRecordSortingIteratorFactory;

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
        final ProgressLogger progressLogger = new ProgressLogger(log, 1000000);
        log.info("Reading in records for TAG name sorting");
        final CloseableIterator<SAMRecord> result =
                SamRecordSortingIteratorFactory.create(writerHeader, reader.iterator(), new StringTagComparator(primaryTag), progressLogger);

		log.info("Sorting finished.");
		return (result);
	}
	

	/**
	 * If the file is sorter in query name order, return an iterator over
	 * the file.  Otherwise, sort records in queryname order and return an iterator over those records.
	 * @return An iterator over the records in the file, in queryname order.
	 */
	
	public static CloseableIterator<SAMRecord> getQuerynameSortedRecords(SamReader reader) {
		if (reader.getFileHeader().getSortOrder().equals(SortOrder.queryname)) {
			return reader.iterator();
		}
		log.info("Input SAM/BAM not in queryname order, sorting...");
		SAMSequenceDictionary dict= reader.getFileHeader().getSequenceDictionary();
		List<SAMProgramRecord> programs =reader.getFileHeader().getProgramRecords();
		
		final SAMFileHeader writerHeader = new SAMFileHeader();
        writerHeader.setSortOrder(SAMFileHeader.SortOrder.queryname);
        writerHeader.setSequenceDictionary(dict);
        for (SAMProgramRecord spr : programs) {
        	writerHeader.addProgramRecord(spr);
        }
        log.info("Reading in records for query name sorting");
        final ProgressLogger progressLogger = new ProgressLogger(log, 1000000, "Sorting reads in query name order");
        final CloseableIterator<SAMRecord> result =
                SamRecordSortingIteratorFactory.create(writerHeader, reader.iterator(), new SAMRecordQueryNameComparator(), progressLogger);
        log.info("Sorting finished.");
        return result; 
	}
	
}
