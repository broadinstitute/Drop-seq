package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.multisample;

import org.testng.Assert;
import org.testng.annotations.Test;



public class SamplePairAssignmentForCellTest {
	@Test
	public void testEqualsHash() {
		SamplePairAssignmentForCell s1 = new SamplePairAssignmentForCell("cell", "s1", "s2", 1d, 0.5d, 0.2d, 0.5d, 1, 2, 10, 1,2,3,4,5);
		SamplePairAssignmentForCell s2 = new SamplePairAssignmentForCell("cell", "s1", "s2", 1d, 0.5d, 0.2d, 0.5d, 1, 2, 10, 1,2,3,4,5);
		
		Assert.assertEquals(s1, s2);
		Assert.assertEquals(s1.toString(), s2.toString());
		s1.hashCode();								
	}
}
