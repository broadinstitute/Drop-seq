package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Interval;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.testng.Assert;
import org.testng.annotations.Test;

public class IntervalTagComparatorTest {

	private final File dictFile = new File ("testdata/org/broadinstitute/transcriptome/utils/hg19.dict");
	private final String intervalTag = "ZI";

	@Test(enabled=true)
	/**
	 * For this test, we generate a lot of data, and sort it.
	 * Since the SAMRecord doesn't really need to change, generate one static record and a LOT of intervals to tag the record with.
	 * Size of the SAMRecord doesn't change, but the speed of the comparator should.
	 */
	public void testSpeed () {
		List<SAMRecord> records = createManyIntervalTaggedSAMRecords(10);


	}

	private List<SAMRecord> createManyIntervalTaggedSAMRecords (final int desiredNumRecords) {
		List<SAMRecord> data = new ArrayList<SAMRecord>();

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
		String s="chr3:10";
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

		List<SAMRecord> result = new ArrayList<SAMRecord>();
		result.add(r1);
		result.add(r2);
		return result;
	}






}
