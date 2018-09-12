package org.broadinstitute.dropseqrna.annotation;

import org.testng.Assert;
import org.testng.annotations.Test;

public class RefFlatRecordTest {

	@SuppressWarnings("unlikely-arg-type")
	@Test
	public void testEquals() {
		RefFlatRecord r1 = new RefFlatRecord("foo", "trans1", "1", false, 1, 10, 2, 9);
		RefFlatRecord r2 = new RefFlatRecord("foo", "trans1", "1", false, 1, 10, 2, 9);
		RefFlatRecord r3 = new RefFlatRecord("foo", "trans2", "2", false, 1, 10, 2, 9);
		r1.addExonStart(2);
	  	r1.addExonEnd(5);
	  	r2.addExonStart(2);
	  	r2.addExonEnd(5);
		boolean flag = r1.equals(r2);
		Assert.assertTrue(flag);
		Assert.assertTrue(r1.equals(r1));

		// add an extra exon to break r1==r2.
		r2.addExonStart(7);
		r2.addExonEnd(8);
		Assert.assertFalse(r1.equals(r2));


		Assert.assertFalse(r1.equals(new String ("nope")));
		Assert.assertFalse(r1.equals(null));
		// test breaking each variable
		r1 = new RefFlatRecord("foo", "trans1", "1", false, 1, 10, 2, 9);
		r2 = new RefFlatRecord("bar", "trans1", "1", false, 1, 10, 2, 9);
		Assert.assertFalse(r1.equals(r2));
		r2 = new RefFlatRecord("foo", "trans2", "1", false, 1, 10, 2, 9);
		Assert.assertFalse(r1.equals(r2));
		r2 = new RefFlatRecord("foo", "trans1", "2", false, 1, 10, 2, 9);
		Assert.assertFalse(r1.equals(r2));
		r1 = new RefFlatRecord("foo", "trans1", "1", true, 1, 10, 2, 9);
		Assert.assertFalse(r1.equals(r2));
		r1 = new RefFlatRecord("foo", "trans1", "1", false, 2, 10, 2, 9);
		Assert.assertFalse(r1.equals(r2));
		r1 = new RefFlatRecord("foo", "trans1", "1", false, 1, 11, 2, 9);
		Assert.assertFalse(r1.equals(r2));
		r1 = new RefFlatRecord("foo", "trans1", "1", false, 1, 10, 3, 9);
		Assert.assertFalse(r1.equals(r2));
		r1 = new RefFlatRecord("foo", "trans1", "1", false, 1, 10, 2, 11);
		Assert.assertFalse(r1.equals(r2));

	}



	@Test
	public void testToString() {
	  	RefFlatRecord r1 = new RefFlatRecord("foo", "trans1", "1", false, 1, 10, 2, 9);
	  	r1.addExonStart(2);
	  	r1.addExonEnd(5);
	  	r1.addExonStart(4);
	  	r1.addExonEnd(7);
	  	String s = r1.toString();
	  	String e = "foo\ttrans1\t1\t+\t0\t10\t1\t9\t2\t1,3,\t5,7,";
	  	Assert.assertNotNull(s);
	  	Assert.assertEquals(s, e);
  }
}
