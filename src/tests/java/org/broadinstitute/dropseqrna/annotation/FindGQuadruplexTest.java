package org.broadinstitute.dropseqrna.annotation;

import java.util.List;

import org.junit.Assert;
import org.testng.annotations.Test;

/**
 * Given a sequence, find a list of GQuadruplex intervals.
 * @author nemesh
 *
 */
public class FindGQuadruplexTest {


	@Test
	public void findSimple () {
		String t1 = "GGGCCTGGGGCTGGGCCTGGG";
		List<GQuadruplex> r1 = GQuadruplex.find("t1", t1);
		Assert.assertEquals(1, r1.size());
		GQuadruplex result = r1.get(0);
		Assert.assertEquals("GGGCCTGGGGCTGGGCCTGGG", result.getSequence());

		Assert.assertEquals(1, result.getMatchInterval().getStart());
		Assert.assertEquals(21, result.getMatchInterval().getEnd());

	}

	@Test
	public void findTandem () {
		String tandem="GGGCCTGGGCCTGGGCCTGGGCCTGGGCCTGGGCCTGGGCCTGGG";
		List<GQuadruplex> r1 = GQuadruplex.find("tandem", tandem);

		Assert.assertEquals(2, r1.size());
		GQuadruplex result = r1.get(0);
		Assert.assertEquals("GGGCCTGGGCCTGGGCCTGGG", result.getSequence());
		Assert.assertEquals(1, result.getMatchInterval().getStart());
		Assert.assertEquals(21, result.getMatchInterval().getEnd());

		result = r1.get(1);
		Assert.assertEquals("GGGCCTGGGCCTGGGCCTGGG", result.getSequence());
		Assert.assertEquals(25, result.getMatchInterval().getStart());
		Assert.assertEquals(45, result.getMatchInterval().getEnd());

	}
}
