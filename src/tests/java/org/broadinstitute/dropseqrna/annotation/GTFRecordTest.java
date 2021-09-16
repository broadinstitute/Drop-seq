package org.broadinstitute.dropseqrna.annotation;

import java.util.List;

import htsjdk.samtools.util.StringUtil;
import org.testng.annotations.Test;

import junit.framework.Assert;

public class GTFRecordTest {

	@Test
	public void testAddError() {
		GTFRecord r1 = new GTFRecord("1", 1, 10, false, "g1ID", "g1,Name", "t1Name", "t1ID", "coding", "exon", 1);
		List<String> errors = r1.validate();
		Assert.assertNotNull(errors);
		Assert.assertEquals(errors.size(), 1);
		String expected="Reserved character ',' in gene name [g1,Name]";
		Assert.assertEquals(expected, errors.get(0));
	}

	@Test
	public void testCompareTo() {
		GTFRecord r1 = new GTFRecord("1", 1, 10, false, "g1ID", "g1Name", "t1Name", "t1ID", "coding", "exon", 1);
		GTFRecord r2 = new GTFRecord("1", 1, 10, false, "g1ID", "g1Name", "t1Name", "t1ID", "coding", "exon", 1);
		Assert.assertTrue(r1.compareTo(r2)==0);
		// only compares on the interval.
		r2 = new GTFRecord("1", 1, 10, false, "g2ID", "g1Name", "t1Name", "t1ID", "coding", "exon", 1);
		Assert.assertTrue(r1.compareTo(r2)==0);
		r2 = new GTFRecord("1", 1, 10, false, "g1ID", "g2Name", "t1Name", "t1ID", "coding", "exon", 1);
		Assert.assertTrue(r1.compareTo(r2)==0);
		r2 = new GTFRecord("1", 1, 10, false, "g1ID", "g1Name", "t2Name", "t1ID", "coding", "exon", 1);
		Assert.assertTrue(r1.compareTo(r2)==0);
		r2 = new GTFRecord("1", 1, 10, false, "g1ID", "g1Name", "t1Name", "t2ID", "coding", "exon", 1);
		Assert.assertTrue(r1.compareTo(r2)==0);
		// change up interval.

		r2 = new GTFRecord("2", 1, 10, false, "g1ID", "g1Name", "t1Name", "t1ID", "coding", "exon", 1);
		int val = r1.compareTo(r2);
		Assert.assertTrue(val<0);

		r2 = new GTFRecord("1", 2, 10, false, "g1ID", "g1Name", "t1Name", "t1ID", "coding", "exon", 1);
		val = r1.compareTo(r2);
		Assert.assertTrue(val<0);
		r2 = new GTFRecord("1", 1, 11, false, "g1ID", "g1Name", "t1Name", "t1ID", "coding", "exon", 1);
		val = r1.compareTo(r2);
		Assert.assertTrue(val<0);
		// strand isn't checked as part of this comparison.
		r2 = new GTFRecord("1", 1, 10, true, "g1ID", "g1Name", "t1Name", "t1ID", "coding", "exon", 1);
		val = r1.compareTo(r2);
		Assert.assertTrue(val==0);


	}

	@Test
	public void testEquals() {
		GTFRecord r1 = new GTFRecord("1", 1, 10, false, "g1ID", "g1Name", "t1Name", "t1ID", "coding", "exon", 1);
		GTFRecord r2 = new GTFRecord("1", 1, 10, false, "g1ID", "g1Name", "t1Name", "t1ID", "coding", "exon", 1);
		Assert.assertTrue(r1.equals(r1));
		Assert.assertFalse(r1.equals(null));
		r2 = new GTFRecord("2", 1, 10, false, "g1ID", "g1Name", "t1Name", "t1ID", "coding", "exon", 1);
		Assert.assertFalse(r1.equals(r2));
		r2 = new GTFRecord("1", 2, 10, false, "g1ID", "g1Name", "t1Name", "t1ID", "coding", "exon", 1);
		Assert.assertFalse(r1.equals(r2));
		r2 = new GTFRecord("1", 1, 11, false, "g1ID", "g1Name", "t1Name", "t1ID", "coding", "exon", 1);
		Assert.assertFalse(r1.equals(r2));
		r2 = new GTFRecord("1", 1, 10, false, "g2ID", "g1Name", "t1Name", "t1ID", "coding", "exon", 1);
		Assert.assertFalse(r1.equals(r2));
		r2 = new GTFRecord("1", 1, 10, false, "g1ID", "g2Name", "t1Name", "t1ID", "coding", "exon", 1);
		Assert.assertFalse(r1.equals(r2));
		r2 = new GTFRecord("1", 1, 10, false, "g1ID", "g1Name", "t2Name", "t1ID", "coding", "exon", 1);
		Assert.assertFalse(r1.equals(r2));
		r2 = new GTFRecord("1", 1, 10, false, "g1ID", "g1Name", "t1Name", "t2ID", "coding", "exon", 1);
		Assert.assertFalse(r1.equals(r2));
		r2 = new GTFRecord("1", 1, 10, false, "g1ID", "g1Name", "t1Name", "t1ID", "intronic", "exon", 1);
		Assert.assertFalse(r1.equals(r2));
		r2 = new GTFRecord("1", 1, 10, false, "g1ID", "g1Name", "t1Name", "t1ID", "coding", "gene", 1);
		Assert.assertFalse(r1.equals(r2));


	}

	@Test
	public void testToString() {
		GTFRecord r1 = new GTFRecord("1", 1, 10, false, "g1ID", "g1,Name", "t1Name", "t1ID", "coding", "exon", 1);
		String actual = r1.toString();
		String expected = "1:1-10\t+\t. [g1,Name exon]";
		Assert.assertEquals(expected, actual);

	}

	@Test
	public void testValidate() {
		// Test no errors
		List<String> errors = new GTFRecord("1", 1, 10, false, "g1ID", "g1Name", "t1Name", "t1ID", "coding", "exon", 1).validate();
		Assert.assertNull(errors);
		// Test null names and IDs
		errors = new GTFRecord("1", 1, 10, false, null, null, null, null, "coding", "exon", 1).validate();
		Assert.assertEquals(StringUtil.join("\n", errors),4, errors.size());
		// Test partially invalid but fixed.
		GTFRecord r = new GTFRecord("1", 1, 10, false, null, "g1Name", "t1Name", null, "coding", "exon", 1);
		r.validate(true);
		Assert.assertEquals(r.getGeneName(), r.getGeneID());
		Assert.assertEquals(r.getTranscriptName(), r.getTranscriptID());
		r = new GTFRecord("1", 1, 10, false, "g1ID", null, null, "t1ID", "coding", "exon", 1);
		r.validate(true);
		Assert.assertEquals(r.getGeneName(), r.getGeneID());
		Assert.assertEquals(r.getTranscriptName(), r.getTranscriptID());

	}

}
