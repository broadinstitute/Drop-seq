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
package org.broadinstitute.dropseqrna.barnyard.digitalexpression;

import org.broadinstitute.dropseqrna.barnyard.digitalexpression.BarcodeSimulator;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.BarcodeSimulator.BarcodeSimulatorResult;
import org.broadinstitute.dropseqrna.utils.editdistance.LevenshteinDistance;
import org.broadinstitute.dropseqrna.utils.editdistance.LevenshteinDistanceResult;
import org.testng.Assert;
import org.testng.Reporter;
import org.testng.annotations.Test;



public class BarcodeSimulatorTest {

	@Test(enabled=true, groups={"dropseq", "transcriptome"})
	public void testSubstition() {
		BarcodeSimulator s = new BarcodeSimulator(12);
		String b1 = s.getRandomBarcode();
		Assert.assertEquals(b1.length(), 12);
		
		String b2=s.addSubstitutions(b1, 2);
		Assert.assertEquals(b1.length(), 12);		
		Assert.assertTrue(!b1.equals(b2));
		
		LevenshteinDistanceResult r =  LevenshteinDistance.computeLevenshteinDistanceResult(b1,b2);
		int ed = r.getEditDistance();
		Assert.assertEquals(ed, 2);
		Assert.assertTrue(!b1.equals(b2));
		
		// for high edit distance barcodes, you can compute the edit distance solution as an insertion or deletion instead of a train of substitutions.
		for (int i=0; i<b1.length(); i++) {
			b2=s.addSubstitutions(b1, i);
			r =  LevenshteinDistance.computeLevenshteinDistanceResult(b1,b2);
			ed = r.getEditDistance();
			String [] ops = r.getOperations();
			// only validate if all operations were substitions.  You can make enough substitutions that an indel inadvertently created has a closer distance.
			int count=0;
			for (String z: ops) {
				if (z.equals("S") || z.equals("M") ) {
					count++;
				}
			}
			// if all operations are "S" or "M" then test the edit distance.
			if (count==b2.length()) {
				Assert.assertEquals(ed, i);
			}
		}
		
	}
	
	@Test(enabled=true, groups={"dropseq", "transcriptome"})
	public void testInsertion() {
		char [] bases = {'A', 'C', 'G', 'T'};
		BarcodeSimulator s = new BarcodeSimulator(12, bases);
		String b1 = s.getRandomBarcode();
		Assert.assertEquals(b1.length(), 12);
		
		BarcodeSimulatorResult b2=s.addInsertion(b1, 1);
		Reporter.log("ED=1 insertion test B1:" + b1 + " B2:" + b2, true);
		Assert.assertEquals(b1.length(), 12);	
		Assert.assertTrue(!b1.equals(b2));
		
		// a single base indel at the start or end of the barcode looks just like a single base substitution.
		// it can also result in no change, say if you have a string like "AAA", and you insert an "A".
		int ed =  LevenshteinDistance.getIndelSlidingWindowEditDistance(b1, b2.modifiedBarcode);
		Assert.assertEquals(ed, 1);
		
		b2=s.addInsertion(b1, 2);
		Reporter.log("ED=2 insertion test  B1:" + b1 + " B2:" + b2, true);
		Assert.assertEquals(b1.length(), 12);	
		Assert.assertTrue(!b1.equals(b2));

		ed =  LevenshteinDistance.getIndelSlidingWindowEditDistance(b1, b2.modifiedBarcode);
		Assert.assertTrue(ed<=2);		
		
	}
	
	@Test(enabled=true, groups={"dropseq", "transcriptome"})
	public void testInsertionAtPosition() {
		int len=12;
		char [] bases = {'A', 'C', 'G', 'T'};
		BarcodeSimulator s = new BarcodeSimulator(len, bases);
		String b1 = s.getRandomBarcode();
		Assert.assertEquals(b1.length(), len);
		
		
		for (int insertionLength=0; insertionLength<len; insertionLength++) {
			for (int position=0; position<len; position++) {
				BarcodeSimulatorResult b2=s.addInsertion(b1, insertionLength, position);
				Assert.assertEquals(b2.modifiedBarcode.length(), len);
			}
			
		}		 
	}
	
	@Test(enabled=true, groups={"dropseq", "transcriptome"})
	public void testSweepInsertions() {
		int len=12;
		char [] bases = {'A', 'C', 'G', 'T'};
		BarcodeSimulator s = new BarcodeSimulator(len, bases);
		String b1 = s.getRandomBarcode();
		Assert.assertEquals(b1.length(), len);
		
		for (int insertionLength=0; insertionLength<=3; insertionLength++) {
			for (int position=0; position<len; position++) {
				for (int i=0; i<100; i++) {
					BarcodeSimulatorResult b2=s.addInsertion(b1, insertionLength);
					Assert.assertEquals(b2.modifiedBarcode.length(), len);
					int ed =  LevenshteinDistance.getIndelSlidingWindowEditDistance(b1, b2.modifiedBarcode);
					if (ed>insertionLength) {
						Reporter.log("insertion length=" + Integer.toString(insertionLength) + " " + b2.toString() + " ed found=" + Integer.toString(ed), true);
					}
					Assert.assertTrue(ed<=insertionLength);	
				}
				
			}
			
		}		 
	}
	
	@Test(enabled=true, groups={"dropseq", "transcriptome"})
	public void testSweepInsertionsPlusSubstitutions() {
		int len=12;
		int maxNumSubs=4;
		int maxNumInsertions=4; 
		int numIterations=1000;
		
		char [] bases = {'A', 'C', 'G', 'T'};
		BarcodeSimulator s = new BarcodeSimulator(len, bases);
		String b1 = s.getRandomBarcode();
		Assert.assertEquals(b1.length(), len);
		
		for (int numSubs=0; numSubs<=maxNumSubs; numSubs++) {
			for (int insertionLength=0; insertionLength<=maxNumInsertions; insertionLength++) {
					for (int i=numIterations; i<1; i++) {
						
						BarcodeSimulatorResult b2=s.addInsertion(b1, insertionLength);
						String finalBC = s.addSubstitutions(b2.modifiedBarcode, numSubs);
						Assert.assertEquals(b2.modifiedBarcode.length(), len);
						int ed =  LevenshteinDistance.getIndelSlidingWindowEditDistance(b1, finalBC);
						if (ed>insertionLength+numSubs) {
							Reporter.log("insertion length=" + Integer.toString(insertionLength) + " " + b2.toString() + " ed found=" + Integer.toString(ed), true);
						}
						Assert.assertTrue(ed<=insertionLength+numSubs);	
					}		
			}	
		}
				 
	}
	
}
