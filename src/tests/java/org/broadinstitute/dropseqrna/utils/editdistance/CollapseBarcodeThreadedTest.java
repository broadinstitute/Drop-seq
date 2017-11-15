/*
 * MIT License
 *
 * Copyright 2017 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package org.broadinstitute.dropseqrna.utils.editdistance;

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
