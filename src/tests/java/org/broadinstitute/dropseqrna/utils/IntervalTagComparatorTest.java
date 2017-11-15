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

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import org.broadinstitute.dropseqrna.utils.readiterators.SamRecordSortingIteratorFactory;
import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.Interval;

public class IntervalTagComparatorTest {

	private final File dictFile = new File ("testdata/org/broadinstitute/transcriptome/utils/hg19.dict");
	private final String intervalTag = "ZI";

	@Test(enabled=false)
	/**
	 * For this test, we generate a lot of data, and sort it.
	 * Since the SAMRecord doesn't really need to change, generate one static record and a LOT of intervals to tag the record with.
	 * Size of the SAMRecord doesn't change, but the speed of the comparator should.
	 */
	public void testSpeed () {
		int numRecords = 30000000;
		int readLength=60;
		SamReader inputSam = SamReaderFactory.makeDefault().open(dictFile);
		SAMSequenceDictionary dict= inputSam.getFileHeader().getSequenceDictionary();
		dict = filterSD(dict);

		long start = System.currentTimeMillis();

		GenerateIntervalTaggedSAMRecords gen = new GenerateIntervalTaggedSAMRecords(dict, intervalTag, readLength, numRecords);
		IntervalTagComparator comparator = new IntervalTagComparator(this.intervalTag);

		final CloseableIterator<SAMRecord> sortingIterator =
	            SamRecordSortingIteratorFactory.create(inputSam.getFileHeader(), gen, comparator, null);

		while(sortingIterator.hasNext()) {
			SAMRecord r = sortingIterator.next();
			Assert.assertNotNull(r);
		}
		long elapsed = System.currentTimeMillis() - start;
		System.out.println("elapsed time = " + elapsed + " ms");
		System.out.println("elapsed time = " + Math.round(elapsed/1000) + " seconds");
	}

	/**
	 * Make this a little more like actual sequence data, where reads mapped to GL are practically non existent.
	 * @param dict
	 * @return
	 */
	private SAMSequenceDictionary filterSD (final SAMSequenceDictionary dict) {
		SAMSequenceDictionary result = new SAMSequenceDictionary();
		for (SAMSequenceRecord r: dict.getSequences())
			if (r.getSequenceLength()>10000000)
				result.addSequence(r);
		return result;
	}


	/**
	 * Replaces reading SAMRecords from a BAM file.  All SAM records are the same, but have a different interval tag on them.
	 * @author nemesh
	 *
	 */
	private class GenerateIntervalTaggedSAMRecords implements Iterator<SAMRecord>{
		private final SAMRecord samRecordTemplate;
		private final List<SAMSequenceRecord> recs;
		private final Random randomGenerator;
		private final String intervalTag;
		private final int readLength;
		private final int maxNumRecords;
		private int currentNumRecords;
		private final int numContigs;

		public GenerateIntervalTaggedSAMRecords (final SAMSequenceDictionary dict, final String intervalTag, final int readLength, final int maxNumRecords) {
			this.intervalTag=intervalTag;
			this.readLength=readLength;
			this.maxNumRecords=maxNumRecords;
			this.recs = dict.getSequences();
			this.numContigs=recs.size();
			randomGenerator = new Random();
			this.currentNumRecords=0;


			final SAMRecordSetBuilder ret = new SAMRecordSetBuilder(false, SAMFileHeader.SortOrder.coordinate);
			ret.addUnmappedFragment("foo");
			this.samRecordTemplate= ret.getRecords().iterator().next();

		}

		@Override
		public SAMRecord next() {
			SAMRecord result = null;
			SAMSequenceRecord r = recs.get(randomGenerator.nextInt(this.numContigs));
			String chr = r.getSequenceName();
			int seqLen = r.getSequenceLength();
			int s = randomGenerator.nextInt(seqLen);
			int e = s+this.readLength-1;
			Interval interval = new Interval (chr, s,e);
			try {
				result = (SAMRecord) samRecordTemplate.clone();
				String tag = IntervalTagComparator.toString(interval);
				// I realize that using encoding the full interval can be a bit heavy handed.
				result.setAttribute(this.intervalTag, tag);

			} catch (CloneNotSupportedException e1) {
				// this should never happen, sigh.
			}
			this.currentNumRecords++;
			return result;
		}

		@Override
		public void remove() {
		}

		@Override
		public boolean hasNext() {
			return (this.currentNumRecords<this.maxNumRecords);
		}

	}

	private List<SAMRecord> createManyIntervalTaggedSAMRecords (final int desiredNumRecords) {
		List<SAMRecord> data = new ArrayList<>();

		SamReader inputSam = SamReaderFactory.makeDefault().open(this.dictFile);
		SAMRecord samRecordTemplate = new SAMRecord (inputSam.getFileHeader());

		SAMSequenceDictionary dict= inputSam.getFileHeader().getSequenceDictionary();
		List<SAMSequenceRecord> recs = dict.getSequences();
		int numRecs = recs.size();

		Random randomGenerator = new Random();
		for (int i=0; i<desiredNumRecords; i++) {
			SAMSequenceRecord r = recs.get(randomGenerator.nextInt(numRecs+1));
			String chr = r.getSequenceName();
			int seqLen = r.getSequenceLength();
			int s1 = randomGenerator.nextInt(seqLen);
			int s2 = randomGenerator.nextInt(seqLen);
			int s = Math.min(s1, s2);
			int e = Math.max(s1, s2);
			Interval interval = new Interval (chr, s1,s2);
			try {
				SAMRecord r1 = (SAMRecord) samRecordTemplate.clone();
				// I realize that using encoding the full interval can be a bit heavy handed.
				r1.setAttribute(this.intervalTag, interval.toString());
				data.add(r1);
			} catch (CloneNotSupportedException e1) {
				// this should never happen, sigh.
			}
		}
		return data;
	}


	@Test
	public void testParseIntervalWithContigContainingDelimiter() {
		Interval i = new Interval ("HLA-A*02:43N", 1205, 1205, false, "foo");
		String s = IntervalTagComparator.toString(i);
		Interval f = IntervalTagComparator.fromString(s);
		Assert.assertEquals(f, i);
	}

	@Test
	public void testParseInterval1() {

		Interval i = new Interval ("chr1", 1, 10, false, "foo");
		String s = IntervalTagComparator.toString(i);
		Interval f = IntervalTagComparator.fromString(s);
		Assert.assertEquals(f, i);
	}

	@Test
	public void testParseInterval2() {
		Interval i = new Interval ("chr1", 1, 10);
		String s = IntervalTagComparator.toString(i);
		Interval f = IntervalTagComparator.fromString(s);
		Assert.assertEquals(f, i);
	}

	@Test
	public void testParseInterval3() {
		// a slightly more ghetto string that's outside the classic interval spec.
		Interval i = new Interval ("chr3", 10, 10);
		String s="chr3|10";
		Interval f = IntervalTagComparator.fromString(s);
		Assert.assertEquals(f, i);
	}


	/**
	 * Tests ordering of a few SAMRecords without use of a dictionary file.
	 * This means that chromosomes are sorted in natural string ordering.
	 * The difference is illustrated below, assuming that the sequence dictionary file has chromosomes
	 * in the order 1,2,3,4,5,6,7,8,9,10,11....so chromosome 2 has index 2, chromosome 10 has index 10.
	 * > sort (c("2", "10"), decreasing=F)
	 * [1] "10" "2"
	 * > sort (c(2, 10), decreasing=F)
	 * [1]  2 10
	 */
	@Test
	public void orderRecordsTestWithDictFile () {
		SamReader inputSam = SamReaderFactory.makeDefault().open(this.dictFile);

		List<SAMRecord> recs = getRecords();

		Integer a = 1;  // these start at 0
		Integer b = 9;
		int exp = a.compareTo(b);

		// when we sort without sequence dictionary, we expect the chr10 tagged read to show up first.
		IntervalTagComparator c = new IntervalTagComparator(this.intervalTag, inputSam.getFileHeader().getSequenceDictionary());
		int comp = c.compare(recs.get(0), recs.get(1));
		// assert sorting same order.
		Assert.assertEquals(exp>0,  comp>0);
	}

	@Test
	public void orderRecordsTestNoDictFile () {
		List<SAMRecord> recs = getRecords();

		String a = "2";
		String b = "10";
		int exp = a.compareTo(b);

		// when we sort without sequence dictionary, we expect the chr10 tagged read to show up first.
		IntervalTagComparator c = new IntervalTagComparator(this.intervalTag);
		int comp = c.compare(recs.get(0), recs.get(1));
		Assert.assertEquals(exp,  comp);
	}

	@Test
	public void testSortIntervalsWithNames () {
		SAMRecordSetBuilder builder  = new SAMRecordSetBuilder();

		// chromosome 2 and chromosome 10
		SAMRecord r1 = builder.addFrag("read1", 1, 1, false);
		Interval i1 = new Interval("2", 1, 10, true, "foo");
		r1.setAttribute(this.intervalTag, IntervalTagComparator.toString(i1));

		SAMRecord r2 = builder.addFrag("read1", 1, 1, false);
		Interval i2 = new Interval("2", 1, 10, true, "bar");
		r2.setAttribute(this.intervalTag, IntervalTagComparator.toString(i2));

		// bar comes before foo.
		int expected = "foo".compareTo("bar");
		IntervalTagComparator c = new IntervalTagComparator(this.intervalTag);
		int comp = c.compare(r1, r2);
		Assert.assertEquals(expected,  comp);
	}



	private List<SAMRecord> getRecords () {
		SAMRecordSetBuilder builder  = new SAMRecordSetBuilder();

		// chromosome 2 and chromosome 10
		SAMRecord r1 = builder.addFrag("read1", 1, 1, false);
		Interval i1 = new Interval("2", 1, 10, true, null);
		r1.setAttribute(this.intervalTag, IntervalTagComparator.toString(i1));

		SAMRecord r2 = builder.addFrag("read1", 1, 1, false);
		Interval i2 = new Interval("10", 1, 10, true, null);
		r2.setAttribute(this.intervalTag, IntervalTagComparator.toString(i2));

		List<SAMRecord> result = new ArrayList<>();
		result.add(r1);
		result.add(r2);
		return result;
	}






}
