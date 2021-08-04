package org.broadinstitute.dropseqrna.sbarro;

import org.testng.annotations.Test;
import org.testng.Assert;

public class BipartiteRabiesVirusCollapseResultTest {

	/**
	 * These are purely for code coverage.
	 */
	@Test
	public void equalsTest() {
		BipartiteRabiesVirusCollapseResult r1 = new BipartiteRabiesVirusCollapseResult("TEST1", "TEST2", 0, 0);
		BipartiteRabiesVirusCollapseResult r2 = new BipartiteRabiesVirusCollapseResult("TEST1", "TEST2", 0, 0, 0.9d);
		Assert.assertEquals(r1, r2);
	}

	@Test
	public void hashCodeTest() {
		BipartiteRabiesVirusCollapseResult r1 = new BipartiteRabiesVirusCollapseResult("TEST1", "TEST2", 0, 0);
		BipartiteRabiesVirusCollapseResult r2 = new BipartiteRabiesVirusCollapseResult("TEST1", "TEST2", 0, 0, 0.9d);
		Assert.assertEquals(r1.hashCode(), r2.hashCode());
	}
}
