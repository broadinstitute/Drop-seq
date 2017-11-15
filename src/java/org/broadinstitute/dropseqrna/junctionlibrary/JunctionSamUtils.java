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
package org.broadinstitute.dropseqrna.junctionlibrary;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Log;

public class JunctionSamUtils {
	private final Log log = Log.getInstance(JunctionSamUtils.class);
	
	private static JunctionSamUtils ourInstance = new JunctionSamUtils();
	
	public static final int MAX_RECORDS_IN_RAM = 500000;
	
	private JunctionSamUtils () {
	}

	/**
	 * Uses singleton pattern so you don't instantiate a million copies with loggers.
	 * @return
	 */
	public static JunctionSamUtils getInstance() {
		return ourInstance;
	}
	
	/**
     * Strips read name of extraneous read number extensions
     */
    public String cleanReadName(String readName) {
        if (readName.endsWith("/1") || readName.endsWith("/2")) {
            readName = readName.substring(0, readName.length()-2);
        }
        return readName;
    }
    
	/**
	 * Test that the records are properly paired.
	 * @param r1 Should be the first read of the pair
	 * @param r2 Should be the second read of the pair
	 */
	public boolean testPairedRead (SAMRecord r1, SAMRecord r2) {
		if (r1==null || r2==null) {
			return (false);
		}
		
		if (r1.getSecondOfPairFlag() || r2.getFirstOfPairFlag()) {
			log.warn(r1.getReadName()+ " passed as 2nd read of pair, should be the first");
			return false;
		}
		if (r2.getFirstOfPairFlag()) {
			log.warn(r2.getReadName()+ " passed as 1st read of pair, should be the second");
			return false;
		}
		boolean result=testSameReadName(r1, r2);
		if (result==false) {
			log.warn("Read names not the same "+ r1.getReadName()+ " "+ r2.getReadName());
		}
		return result;
	}
	
	/**
	 * Reads share the same name, and one is flagged first of pair, and the other 2nd of pair.
	 * Order of reads is unimportant.
	 * @param r1
	 * @param r2
	 * @return
	 */
	public boolean testReadsArePaired (SAMRecord r1, SAMRecord r2) {
		// test if the pair of reads are flagged first and second read.
		boolean f1=r1.getFirstOfPairFlag() & r2.getSecondOfPairFlag();
		boolean f2=r1.getSecondOfPairFlag() & r2.getFirstOfPairFlag();
		if (f1 | f2) {
			return (testSameReadName(r1,r2));
		}
		return false;
	}
	
	public boolean testSameReadName (SAMRecord r1, SAMRecord r2) {
		if (r1==null || r2==null) return (false);
		String cleanReadNameOne=cleanReadName(r1.getReadName());
		String cleanReadNameTwo=cleanReadName(r2.getReadName());
		if (cleanReadNameOne.equals(cleanReadNameTwo)==false) {
			return false;
		}
		return true;
	}
	
	/**
	 * If the file is sorter in query name order, return an iterator over
	 * the file.  Otherwise, sort records in queryname order and return an iterator over those records.
	 * @param bamFile The BAM/SAM file to read
	 * @return An iterator over the records in the file, in queryname order.
	 */
	/*
	public CloseableIterator<SAMRecord> getQuerynameSortedRecords(SamReader reader) {
		CloseableIterator<SAMRecord> iter = null;
		if (reader.getFileHeader().getSortOrder().equals(SortOrder.queryname)) {
			iter = reader.iterator();
			return iter;
		}
		log.info("Input SAM/BAM not in queryname order, sorting...");
		ProgressLogger p = new ProgressLogger(this.log, 1000000, "Sorting reads in query name order");
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
	*/

	
    
}
