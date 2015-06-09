package org.broadinstitute.dropseqrna.priv.utils.editdistance;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.apache.commons.lang.RandomStringUtils;
import org.broadinstitute.dropseqrna.utils.editdistance.EDUtils;
import org.broadinstitute.dropseqrna.utils.editdistance.LevenshteinDistance;
import org.broadinstitute.dropseqrna.utils.editdistance.LevenshteinDistanceResult;
import org.junit.Assert;
import org.testng.annotations.Test;

public class CollapseBarcodeThreadedTest {

	@Test(enabled=true, groups={"transcriptome", "dropseq"})
	public void test() {
		String baseString="democrat";
		List<String> comparisonStrings = new ArrayList<String>();
		comparisonStrings.add("republican");
		
		LevenshteinDistanceResult r =  LevenshteinDistance.computeLevenshteinDistanceResult(baseString, comparisonStrings.get(0));
		int ed = r.getEditDistance();
		Assert.assertEquals(8, ed);
		
	}
	
	private List<String> getRandomBarcodes (int numBases, int numBarcodes) {
		char [] bases = {'A', 'C', 'G', 'T', 'N'};
		List<String> result = new ArrayList<String>(numBarcodes);
		for (int i=0; i<numBarcodes; i++) {
			result.add(RandomStringUtils.random(numBases, bases));
		}
		return (result);
		
	}
	
	@Test(enabled=false, groups={"transcriptome", "dropseq"})
	public void test2 () {
		int numBases=8;
		int numBarcodes=10000;
		int chunkSize=1000;
		int editDistance=2;
		int threshold=5;
		
		
		String baseString = getRandomBarcodes(numBases, 1).get(0);
		List<String> comparisonStrings = getRandomBarcodes(numBases, numBarcodes);
		long startTime = System.currentTimeMillis();
		Set<String> closeBarcodes = EDUtils.getInstance().getStringsWithinEditDistanceWithIndel(baseString,comparisonStrings, editDistance);
		Assert.assertNotNull(closeBarcodes);
		long endTime = System.currentTimeMillis();
		long duration = endTime - startTime;
		System.out.println("Single threaded original method took [" + duration + "] miliseconds to process [" + comparisonStrings.size()+ "] barcodes");
		
		long startTime2 = System.currentTimeMillis();
		CollapseBarcodeThreaded cbt = new CollapseBarcodeThreaded(chunkSize, null);
		Set<String> closeBarcodes2 = cbt.getStringsWithinEditDistanceWithIndel(baseString, comparisonStrings, editDistance, true);
		Assert.assertNotNull(closeBarcodes2);
		long endTime2 = System.currentTimeMillis();
		long duration2 = endTime2 - startTime2;
		System.out.println("Multithreaded threaded method with blocksize ["+chunkSize+"] took " + duration2 + " miliseconds and found [" + closeBarcodes2.size() + "] nearby barcodes" );
		double speedBoost=(double)duration/(double)duration2;
		
		Assert.assertEquals(closeBarcodes, closeBarcodes2);
		System.out.println("Speed increase [" + Double.toString(speedBoost) + "] number of threads used ["+ cbt.getNumThreads()+"]");
		
	}
	
	@Test(enabled=false, groups={"transcriptome", "dropseq"})
	public void test3 () {
		int numBases=8;
		int numBaseStrings=1000;
		int numBarcodes=10000;
		int chunkSize=1000;
		int editDistance=1;
		
		List<String> baseStrings = getRandomBarcodes(numBases, numBaseStrings);
		List<String> comparisonStrings = getRandomBarcodes(numBases, numBarcodes);
		
		long startTime = System.currentTimeMillis();
		CollapseBarcodeThreaded cbt = new CollapseBarcodeThreaded(chunkSize, null);
		for (String baseString: baseStrings) {
			Set<String> closeBarcodes = cbt.getStringsWithinEditDistanceWithIndel(baseString, comparisonStrings, editDistance, true);
			Assert.assertNotNull(closeBarcodes);
		}
		
		
		long endTime = System.currentTimeMillis();
		long duration = endTime - startTime;
		System.out.println("Multithreaded threaded method with blocksize ["+chunkSize+"] over ["+ baseStrings.size()+"] strings against [" + comparisonStrings.size() + "] comparison strings took " + duration + " miliseconds for ["+ numBaseStrings*numBarcodes +"] total comparisons");
		
	}

}
