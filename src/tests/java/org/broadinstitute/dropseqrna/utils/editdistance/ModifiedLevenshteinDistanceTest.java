package org.broadinstitute.dropseqrna.utils.editdistance;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.dropseqrna.utils.editdistance.LevenshteinDistance;
import org.broadinstitute.dropseqrna.utils.editdistance.LevenshteinDistanceResult;
import org.junit.Assert;
import org.testng.annotations.Test;

public class ModifiedLevenshteinDistanceTest {
	
	@Test(enabled=true, groups={"transcriptome", "dropseq"})
	public void websiteTest () {
		// see http://stackoverflow.com/questions/5849139/levenshtein-distance-inferring-the-edit-operations-from-the-matrix
		
		LevenshteinDistanceResult r =  LevenshteinDistance.computeLevenshteinDistanceResult("democrat", "republican");
		int ed = r.getEditDistance();
		Assert.assertEquals(8, ed);
		String [] operations = r.getOperations();
		String result = StringUtils.join(operations);
		Assert.assertEquals(result, "SMIISSSSMS");
		
	}
	
	@Test(enabled=true, groups={"transcriptome", "dropseq"})
	public void simpleTests1 () {
		LevenshteinDistanceResult r =  LevenshteinDistance.computeLevenshteinDistanceResult("AAA", "AAA");
		int ed = r.getEditDistance();
		Assert.assertEquals(0, ed);
		
		r =  LevenshteinDistance.computeLevenshteinDistanceResult("AAA", "ATA");
		ed = r.getEditDistance();
		Assert.assertEquals(1, ed);
		
		r =  LevenshteinDistance.computeLevenshteinDistanceResult("AAA", "ATA", 1,1,2);
		ed = r.getEditDistance();
		Assert.assertEquals(2, ed);
		
		// test deletion
		r =  LevenshteinDistance.computeLevenshteinDistanceResult("AATA", "AAA");
		ed = r.getEditDistance();
		Assert.assertEquals(1, ed);
		
		// test insertion
		r =  LevenshteinDistance.computeLevenshteinDistanceResult("AAA", "AATA");
		ed = r.getEditDistance();
		Assert.assertEquals(1, ed);
		
		// test deletion - modify deletion weight
		r =  LevenshteinDistance.computeLevenshteinDistanceResult("AATA", "AAA", 2,1,1);
		ed = r.getEditDistance();
		Assert.assertEquals(2, ed);
				
		// test insertion - modify insertion weight
		r =  LevenshteinDistance.computeLevenshteinDistanceResult("AAA", "AATA", 1,2,1);
		ed = r.getEditDistance();
		Assert.assertEquals(2, ed);
		
		
	}
	
	@Test(enabled=true, groups={"transcriptome", "dropseq"})
	public void testEditDistanceOperations () {
		LevenshteinDistanceResult r =  LevenshteinDistance.computeLevenshteinDistanceResult("AAA", "AATA");
		int ed = r.getEditDistance();
		Assert.assertEquals(1, ed);
		String [] operations = r.getOperations();
		String result = StringUtils.join(operations);
		Assert.assertEquals("MMIM", result);
		
	}
	
	
	@Test(enabled=true, groups={"transcriptome", "dropseq"})
	public void testEditDistanceOperationsFull () {
		
		testEditDistance("AGGCTAGGTGTT","AGGCTAAGGTGT", 1);
		testEditDistance("AGGCTAGGTGTT","AGCTAGGTGTTT",1);
		testEditDistance("AGGCTAGGTGTT","AGCTAGCTGTTT",2);
		
		// 1 indel, 2 subs
		testEditDistance("AGGCTAGGTGTT","ACCTAGCTGTTT",3);
		
		// 1 indel, 0 subs
		testEditDistance("AAAGAGGACGG","AAGAGGACGCG",1);
				
		
	}
	
	@Test(enabled=true, groups={"transcriptome", "dropseq"})
	public void testEditDistanceOperationsFull2 () {
		
		// 2 bp insertion
		testEditDistance("GCT", "TAG", 2);
		
		// 2 bp substitution
		testEditDistance("AATG","ZZTG",2);
		
		// 2 bp insertion of a "GC"
		testEditDistance("TCCCAGCAATAG","TCCCAGCAAGCT",2);
		
		// has an indel of length 1 and a substitution of the first base.
		testEditDistance("TATTTTACGCAG","AATTTTAACGCA",2);
		
		// 2 bp insertion, 1 substiution
		testEditDistance("ACCCAGCAATAG","TCCCAGCAAGCT",3);
		
		// go bigger.
		
		// 4 bp insertion, 2 substiution.  Turns out this is ED=4, 'cause long runs of a base "lose" info. 
		testEditDistance("TCCCAGCAAGCT","TCCAAAACATTA",4);
		
		// 3 bp insertion, 3bp substitution.
		testEditDistance("AAACCATTTGGG","GAGCCCCCCTTT",6);
		
		
		
	}
	
	private void testEditDistance (String a, String b, int expectedED) {
		LevenshteinDistanceResult r = LevenshteinDistance.computeLevenshteinDistanceResult(a,b, 1,1,2);
		
		int ed2 = r.getEditDistanceIndelCorrected();
		Assert.assertEquals(expectedED, ed2);
		
	}
	
	
}
