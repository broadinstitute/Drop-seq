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

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import junit.framework.Assert;

import org.broadinstitute.dropseqrna.utils.FilterBAMByTag;
import org.testng.annotations.Test;

public class FilterBAMByTagTest {

	@Test(enabled = true, groups = { "dropseq", "transcriptome" })
	public void filterReadTest() {
		SAMRecord readHasAttribute = new SAMRecord(null);
		String tag = "XT";
		readHasAttribute.setAttribute(tag, "1");

		Set<String> values = new HashSet<String>();
		values.add("1");

		SAMRecord readNoAttribute = new SAMRecord(null);

		FilterBAMByTag t = new FilterBAMByTag();
		// read has attribute, accept any value, want to retain read.
		boolean flag1 = t.filterRead(readHasAttribute, tag, null, true);
		Assert.assertFalse(flag1);

		// read has attribute, accept any value, want to filter read.
		boolean flag2 = t.filterRead(readHasAttribute, tag, null, false);
		Assert.assertTrue(flag2);

		// read has attribute, accept certain value, want to retain read.
		boolean flag3 = t.filterRead(readHasAttribute, tag, values, true);
		Assert.assertFalse(flag3);

		// read has attribute, accept certain value, want to filter read.
		boolean flag4 = t.filterRead(readHasAttribute, tag, values, false);
		Assert.assertTrue(flag4);

		// read does not have attribute, accept any value, want to retain read.
		boolean flag5 = t.filterRead(readNoAttribute, tag, null, true);
		Assert.assertTrue(flag5);

		// read does not have attribute, accept any value, want to filter read.
		boolean flag6 = t.filterRead(readNoAttribute, tag, null, false);
		Assert.assertFalse(flag6);

		// read does not have attribute, accept certain value, want to retain read.
		boolean flag7 = t.filterRead(readNoAttribute, tag, values, true);
		Assert.assertTrue(flag7);

		// read does not have attribute, accept certain value, want to filter read.
		boolean flag8 = t.filterRead(readNoAttribute, tag, values, false);
		Assert.assertFalse(flag8);

	}

	/**
	 * Returns a paired read, first of pair in the first position of the list, 2nd of pair in the 2nd position.
	 * @return
	 */
	private List<SAMRecord> getPairedRead () {
		List<SAMRecord> result = new ArrayList<SAMRecord> ();

		SAMRecordSetBuilder builder = new SAMRecordSetBuilder();
		builder.addUnmappedPair("test");
		Collection<SAMRecord> recs = builder.getRecords();

		for (SAMRecord r: recs) {
			if (r.getFirstOfPairFlag()) result.add(0, r);
			if (r.getSecondOfPairFlag()) result.add(1, r);
		}
		return (result);

	}

	@Test(enabled = true, groups = { "dropseq", "transcriptome" })
	public void filterByReadNumberTest() {
		FilterBAMByTag t = new FilterBAMByTag();

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

}
