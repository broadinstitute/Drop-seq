package org.broadinstitute.dropseqrna.annotation;

import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.util.Interval;

public class GQuadruplexTest {

	@Test
	public void f() {
		// handy examples http://bsbe.iiti.ac.in/bsbe/ipdb/g4dna.php
		// I trimmed a set of starting T's off the 5' end.
		String seq1 = "GGGGACTTTCCGGGAGGCGTGGGGGTTTTTGGGGG";

		// not a G-quadraplex.  Random sequence from the start of chromosome 1, hg19.
		String seq4 = "GGACGCATTTAAAGCAGTGTGTAAAGAGACATTTATAGCACTAAATGCCCACAAGAGACCTCTGCCTGAGAACGTGGGTTTCAGCCTAAGAGTTGTAATA";

		List<GQuadruplex> r1 = GQuadruplex.find("valid1", seq1);
		List<GQuadruplex> r4 = GQuadruplex.find("invalid1", seq4);
		Assert.assertNotNull(r1);
		Assert.assertTrue(r1.size()==1);
		GQuadruplex t1= r1.get(0);  // only 1 is found.

		Assert.assertEquals(35, t1.getMatchInterval().length());
		Assert.assertEquals(t1.getMatchInterval(), new Interval ("valid1", 1,35));
		//G1: GGGG
		//L1: ACTTTCC
		//G2: GGG
		//L2: AGGCGTG
		//G3: GGGG
		//L3: TTTTTGG
		//G4: GGG
		Assert.assertEquals(t1.getG1(), "GGGG");
		Assert.assertEquals(t1.getG2(), "GGG");
		Assert.assertEquals(t1.getG3(), "GGGG");
		Assert.assertEquals(t1.getG4(), "GGG");
		Assert.assertEquals(t1.getL1(), "ACTTTCC");
		Assert.assertEquals(t1.getL2(), "AGGCGTG");
		Assert.assertEquals(t1.getL3(), "TTTTTGG");
		t1.toString();

	}
}
