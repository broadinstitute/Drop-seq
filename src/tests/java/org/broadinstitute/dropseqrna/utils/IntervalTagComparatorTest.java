package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.Interval;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

public class IntervalTagComparatorTest {

	private final File dictFile = new File ("testdata/org/broadinstitute/transcriptome/utils/hg19.dict");
	private final String intervalTag = "ZI";

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
