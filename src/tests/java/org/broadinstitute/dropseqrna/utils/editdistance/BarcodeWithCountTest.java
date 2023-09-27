package org.broadinstitute.dropseqrna.utils.editdistance;

import org.testng.annotations.Test;

import org.testng.Assert;

public class BarcodeWithCountTest {
	@SuppressWarnings("unlikely-arg-type")
	@Test
	public void testBarcodeWithCountTest() {

		BarcodeWithCount b1 = new BarcodeWithCount("AAAA", 5);
		BarcodeWithCount b2 = new BarcodeWithCount("AAAA", 2);
		BarcodeWithCount b3 = new BarcodeWithCount("GGGG", 3);
		BarcodeWithCount b4 = new BarcodeWithCount("GGGG", 7);

		Assert.assertEquals("AAAA", b1.getBarcode());
		Assert.assertSame(5, b1.getCount());

		// equals uses barcode name and count.
		Assert.assertTrue(b1.equals(b1));
		Assert.assertFalse(b1.equals(b2));
		Assert.assertFalse(b1.equals(b3));

		// other dumb tests for coverage
		Assert.assertFalse(b1.equals(null));
		Assert.assertFalse(b1.equals(new String ("Foo")));
		Assert.assertNotNull(b1.hashCode());

		// compares by count alone.
		BarcodeWithCount.CountComparator c = new BarcodeWithCount.CountComparator();
		int r = c.compare(b1, b2);
		Assert.assertTrue(r>0);




	}
}
