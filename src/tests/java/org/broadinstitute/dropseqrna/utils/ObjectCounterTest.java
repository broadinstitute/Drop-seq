package org.broadinstitute.dropseqrna.utils;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.testng.annotations.Test;

import junit.framework.Assert;

public class ObjectCounterTest {
	@Test
	public void test() {
		ObjectCounter<String> o = new ObjectCounter<>();
		o.increment("FOO");
		o.incrementByCount("BAR", 4);

		ObjectCounter<String> o2 = new ObjectCounter<>();
		o2.increment("ZOO");
		o2.incrementByCount("ZAR", 4);

		ObjectCounter<String> both = new ObjectCounter<>(o);
		both.increment(o2);

		both.decrement("BAR");  // now 3.
		both.decrementByCount("ZAR", 2); // now 2.

		Set<String> expectedKeys = new HashSet<>(Arrays.asList("FOO", "BAR", "ZOO", "ZAR"));
		Assert.assertEquals(expectedKeys, both.getKeys());

		List<String> expectedOrderedKeys = Arrays.asList("BAR", "ZAR", "FOO", "ZOO");
		List<String> orderedKeys = both.getKeysOrderedByCount(true);
		Assert.assertEquals(expectedOrderedKeys, orderedKeys);

		// notice keys are ordered by natural order when there are ties to counts in FOO and ZOO.
		List<String> expectedOrderedKeys2 = Arrays.asList( "FOO", "ZOO", "ZAR", "BAR");
		List<String> orderedKeys2 = both.getKeysOrderedByCount(false);
		Assert.assertEquals(expectedOrderedKeys2, orderedKeys2);

		Assert.assertSame(4, both.getSize());
		Assert.assertSame(3, both.getCountForKey("BAR"));
		Assert.assertSame(7, both.getTotalCount());
		Assert.assertSame(2, both.getNumberOfSize(1));
		Assert.assertEquals("BAR",  both.getMode()); // bar is the most common, AKA the highest count.
		both.increment("ZOO");
		Assert.assertEquals("FOO",  both.getMin());

		both.filterByMinCount(2);
		Assert.assertSame(3, both.getSize());

		Assert.assertTrue (both.hasKey("ZOO"));
		Assert.assertFalse (both.hasKey("ZOOPPP"));

		both.setCount("ZOOPPP", 8);
		Assert.assertSame(8, both.getCountForKey("ZOOPPP"));
		both.remove("ZOOPPP");
		Assert.assertSame(0, both.getCountForKey("ZOOPPP"));

		String expected = "{BAR=3, ZAR=2, ZOO=2}";
		String actual = both.toString();
		Assert.assertEquals(expected, actual);

		both.clear();
		Assert.assertSame(0, both.getCounts().size());





	}
}
