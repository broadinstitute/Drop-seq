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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import htsjdk.samtools.metrics.MetricsFile;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class FilterBamByTagTest {

	private static final int PAIRED_READS_ACCEPTED = 12;
	private static final int PAIRED_READS_REJECTED = 8;
	private static final int UNPAIRED_READS_ACCEPTED = 6;
	private static final int UNPAIRED_READS_REJECTED = 4;
	private static final File PAIRED_INPUT_FILE=new File ("testdata/org/broadinstitute/dropseq/utils/paired_reads_tagged.bam");
	private static final File UNPAIRED_INPUT_FILE=new File ("testdata/org/broadinstitute/dropseq/utils/unpaired_reads_tagged.bam");
	private static final File PAIRED_INPUT_FILE_FILTERED=new File ("testdata/org/broadinstitute/dropseq/utils/paired_reads_tagged_filtered.bam");
	private static final File UNPAIRED_INPUT_FILE_FILTERED=new File ("testdata/org/broadinstitute/dropseq/utils/unpaired_reads_tagged_filtered.bam");
	private static final File PAIRED_INPUT_CELL_BARCODES=new File ("testdata/org/broadinstitute/dropseq/utils/paired_reads_tagged.cell_barcodes.txt");

	private static final File UNPAIRED_INPUT_FILE_FILTERED_AAAGTAGAGTGG=new File ("testdata/org/broadinstitute/dropseq/utils/unpaired_reads_tagged_filtered_AAAGTAGAGTGG.bam");

	/**
	 * Runs the CLP, in the given mode, and asserts success
	 * @param successThreshold If non-null, passed to PASSING_READ_THRESHOLD CLP argument
	 * @return the CLP object after doWork() called.
	 */
	private FilterBamByTag runClp(final boolean pairedMode, final Double successThreshold) throws IOException {
		FilterBamByTag f = new FilterBamByTag();
		final String prefix;
		if (pairedMode) {
			f.INPUT=PAIRED_INPUT_FILE;
			f.TAG_VALUES_FILE=PAIRED_INPUT_CELL_BARCODES;
			prefix = "paired_input";
		} else {
			f.INPUT=UNPAIRED_INPUT_FILE;
			// For some reason, use the same file as for paired
			f.TAG_VALUES_FILE = PAIRED_INPUT_CELL_BARCODES;
			prefix = "unpaired_input";
		}
		f.TAG="XC";
		f.PAIRED_MODE = pairedMode;
		f.PASSING_READ_THRESHOLD = successThreshold;
		f.OUTPUT=File.createTempFile(prefix, ".bam");
		f.OUTPUT.deleteOnExit();
		f.SUMMARY=File.createTempFile(prefix, ".summary.txt");
		f.SUMMARY.deleteOnExit();
		Assert.assertEquals(f.doWork(), 0);
		return f;
	}

	@Test
	public void testDoWorkPaired () throws IOException {
		final FilterBamByTag f = runClp(true, null);
		//samtools view -c paired_reads_tagged_filtered.bam
		// 12
		// samtools view -c paired_reads_tagged.bam
		// 20

		List<FilteredReadsMetric> metrics = MetricsFile.readBeans(f.SUMMARY);
		Assert.assertEquals(PAIRED_READS_ACCEPTED, metrics.get(0).READS_ACCEPTED);
		Assert.assertEquals(PAIRED_READS_REJECTED, metrics.get(0).READS_REJECTED);
								
		CompareBAMTagValues cbtv = new CompareBAMTagValues();
		cbtv.INPUT_1=PAIRED_INPUT_FILE_FILTERED;
		cbtv.INPUT_2=f.OUTPUT;
		List<String> tags = new ArrayList<>();
		tags.add("XC");
		cbtv.TAGS=tags;
		int r = cbtv.doWork();
		Assert.assertEquals(r, 0);
		
	}

	@Test
	public void testDoWorkUnPaired () throws IOException {
		FilterBamByTag f = runClp(false, null);

		// samtools view -c unpaired_reads_tagged_filtered.bam
		// 6
		// samtools view -c unpaired_reads_tagged.bam 
		// 10
				
		List<FilteredReadsMetric> metrics = MetricsFile.readBeans(f.SUMMARY);
		Assert.assertEquals(UNPAIRED_READS_ACCEPTED, metrics.get(0).READS_ACCEPTED);
		Assert.assertEquals(UNPAIRED_READS_REJECTED, metrics.get(0).READS_REJECTED);

		CompareBAMTagValues cbtv = new CompareBAMTagValues();
		cbtv.INPUT_1=UNPAIRED_INPUT_FILE_FILTERED;
		cbtv.INPUT_2=f.OUTPUT;
		List<String> tags = new ArrayList<>();
		tags.add("XC");
		cbtv.TAGS=tags;
		int r = cbtv.doWork();
		Assert.assertEquals(r, 0);

		// test alternate path without tag values file.
		f.INPUT=UNPAIRED_INPUT_FILE;
		f.OUTPUT=File.createTempFile("unpaired_input_single_cell", ".bam");
		f.TAG="XC";
		f.TAG_VALUE="AAAGTAGAGTGG";
		f.TAG_VALUES_FILE=null;
		f.PAIRED_MODE=false;
		f.OUTPUT.deleteOnExit();
		Assert.assertEquals(f.doWork(), 0);


		cbtv.INPUT_1=UNPAIRED_INPUT_FILE_FILTERED_AAAGTAGAGTGG;
		cbtv.INPUT_2=f.OUTPUT;
		cbtv.TAGS=tags;
		r = cbtv.doWork();
		Assert.assertEquals(r, 0);
		

	}


	@Test
	public void filterReadTest() {
		SAMRecord readHasAttribute = new SAMRecord(null);
		String tag = "XT";
		readHasAttribute.setAttribute(tag, "1");

		Set<String> values = new HashSet<>();
		values.add("1");

		SAMRecord readNoAttribute = new SAMRecord(null);

		FilterBamByTag t = new FilterBamByTag();
		// read has attribute, accept any value, want to retain read.
		boolean flag1 = t.filterRead(readHasAttribute, tag, null, true, null);
		Assert.assertFalse(flag1);

		// read has attribute, accept any value, want to filter read.
		boolean flag2 = t.filterRead(readHasAttribute, tag, null, false, null);
		Assert.assertTrue(flag2);

		// read has attribute, accept certain value, want to retain read.
		boolean flag3 = t.filterRead(readHasAttribute, tag, values, true, null);
		Assert.assertFalse(flag3);

		// read has attribute, accept certain value, want to filter read.
		boolean flag4 = t.filterRead(readHasAttribute, tag, values, false, null);
		Assert.assertTrue(flag4);

		// read does not have attribute, accept any value, want to retain read.
		boolean flag5 = t.filterRead(readNoAttribute, tag, null, true, null);
		Assert.assertTrue(flag5);

		// read does not have attribute, accept any value, want to filter read.
		boolean flag6 = t.filterRead(readNoAttribute, tag, null, false, null);
		Assert.assertFalse(flag6);

		// read does not have attribute, accept certain value, want to retain read.
		boolean flag7 = t.filterRead(readNoAttribute, tag, values, true, null);
		Assert.assertTrue(flag7);

		// read does not have attribute, accept certain value, want to filter read.
		boolean flag8 = t.filterRead(readNoAttribute, tag, values, false, null);
		Assert.assertFalse(flag8);
		
		// test map quality filtering
		
		readHasAttribute.setMappingQuality(10);
		boolean flag9 = t.filterRead(readHasAttribute, tag, null, true, 10);
		Assert.assertFalse(flag9);
		boolean flag10 = t.filterRead(readHasAttribute, tag, null, true, 20);
		Assert.assertTrue(flag10);
		

	}

	/**
	 * @return a paired read, first of pair in the first position of the list, 2nd of pair in the 2nd position.
	 */
	private List<SAMRecord> getPairedRead () {
		List<SAMRecord> result = new ArrayList<> ();

		SAMRecordSetBuilder builder = new SAMRecordSetBuilder();
		builder.addUnmappedPair("test");
		Collection<SAMRecord> recs = builder.getRecords();

		for (SAMRecord r: recs) {
			if (r.getFirstOfPairFlag()) result.add(0, r);
			if (r.getSecondOfPairFlag()) result.add(1, r);
		}
		return (result);

	}

	@Test
	public void filterByReadNumberTest() {
		FilterBamByTag t = new FilterBamByTag();

		// record paired and read is 1st
		List<SAMRecord> recs = getPairedRead ();
		SAMRecord recFirstPaired = recs.get(0);
		SAMRecord recSecondPaired = recs.get(1);

		boolean flag1= t.retainByReadNumber(recFirstPaired, 1);
		boolean flag2= t.retainByReadNumber(recFirstPaired, 2);
		Assert.assertTrue(flag1);
		Assert.assertFalse(flag2);

		// record paired and read is 2st
		recSecondPaired.setProperPairFlag(true);
		recSecondPaired.setSecondOfPairFlag(true);
		flag1= t.retainByReadNumber(recSecondPaired, 1);
		flag2= t.retainByReadNumber(recSecondPaired, 2);
		Assert.assertTrue(flag2);
		Assert.assertFalse(flag1);

		// record unpaired and read is 1st
		SAMRecordSetBuilder builder = new SAMRecordSetBuilder();
		builder.addUnmappedFragment("foo");
		SAMRecord recFirstUnPaired = builder.getRecords().iterator().next();

		flag1= t.retainByReadNumber(recFirstUnPaired, 1);
		flag2= t.retainByReadNumber(recFirstPaired, 2);
		Assert.assertTrue(flag1);
		Assert.assertFalse(flag2);
	}

	@Test
	public void testArgErrors () throws IOException {
		FilterBamByTag f = new FilterBamByTag();
		f.INPUT=PAIRED_INPUT_FILE;
		f.OUTPUT=File.createTempFile("paired_input", ".bam");
		f.PAIRED_MODE=true;
		f.OUTPUT.deleteOnExit();
		Assert.assertSame(1, f.doWork());

	}

}
