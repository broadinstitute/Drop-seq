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

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import junit.framework.Assert;

import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.editdistance.EDUtils;
import org.testng.annotations.Test;

public class MapBarcodesByEditDistanceTest {

	private static File testData = new File ("testdata/org/broadinstitute/transcriptome/utils/editdistance/inEditDistSmall.txt");
	
	@Test(enabled=true)
	public void collapseBarcodesLarge() {
		int numCoreCells=50;
		
		ObjectCounter <String> barcodes = EDUtils.readBarCodeFile(testData);
		MapBarcodesByEditDistance mapper = new  MapBarcodesByEditDistance(true, 4, 0);
		List<String> barcodeStrings = barcodes.getKeysOrderedByCount(true);
		List<String> coreBarcodes = barcodeStrings.subList(0, numCoreCells);
		
		Map<String, List<String>> result = mapper.collapseBarcodes(coreBarcodes, barcodes, false, 1);
		Assert.assertNotNull(result);
		for (String key : result.keySet()) {
			List<String> values = result.get(key);
			checkResult(key, values);
		}
		
		
		
		
	}
	
	private void checkResult (String key, List<String> values) {
		// TTTTTTTTTTTT=[ATTTTTTTTTTT, CTTTTTTTTTTT, GTTTTTTTTTTT, TATTTTTTTTTT, TCTTTTTTTTTT, TGTTTTTTTTTT, TTTTTTTTTTCT, TTTTTTTTTTTC]
		if (key.equals("TTTTTTTTTTTT")) {
			String [] expected = {"ATTTTTTTTTTT", "CTTTTTTTTTTT", "GTTTTTTTTTTT", "TATTTTTTTTTT", "TCTTTTTTTTTT", "TGTTTTTTTTTT", "TTTTTTTTTTCT", "TTTTTTTTTTTC"};
			checkResult(key, values, expected);
		}
		//AATTTTTTTTTT=[AAATTTTTTTTT, AAGTTTTTTTTT, AATTTTTTTTTA, AATTTTTTTTTC, ACTTTTTTTTTT, AGTTTTTTTTTT, CATTTTTTTTTT, GATTTTTTTTTT]
		if (key.equals("AATTTTTTTTTT")) {
			String [] expected = {"AAATTTTTTTTT", "AAGTTTTTTTTT", "AATTTTTTTTTA", "AATTTTTTTTTC", "ACTTTTTTTTTT", "AGTTTTTTTTTT", "CATTTTTTTTTT", "GATTTTTTTTTT"};
			checkResult(key, values, expected);
		}
		// AAGGGGGGGGGG=[AAAGGGGGGGGG, ATGGGGGGGGGG, TAGGGGGGGGGG]
		if (key.equals("AAGGGGGGGGGG")) {
			String [] expected = {"AAAGGGGGGGGG", "ATGGGGGGGGGG", "TAGGGGGGGGGG"};
			checkResult(key, values, expected);
		}
		// GGGGCGGGGGGA=[GGGGAGGGGGGA, GGGGCGGGGGCA, GGGGCGGGGGGC, GTGGCGGGGGGA]
		if (key.equals("GGGGCGGGGGGA")) {
			String [] expected = {"GGGGAGGGGGGA", "GGGGCGGGGGCA", "GGGGCGGGGGGC", "GTGGCGGGGGGA"};
			checkResult(key, values, expected);
		}
		// GGGGGGGGGGTT=[AGGGGGGGGGTT, CGGGGGGGGGTT, GGGGGGGGGGAT, GGGGGGGGGTTT, TGGGGGGGGGTT] 
		if (key.equals("GGGGGGGGGGTT")) {
			String [] expected = {"AGGGGGGGGGTT", "CGGGGGGGGGTT", "GGGGGGGGGGAT", "GGGGGGGGGTTT", "TGGGGGGGGGTT"};
			checkResult(key, values, expected);
		}
		// CGGCGGGGGGGG=[AGGCGGGGGGGG, CGGCGGGGGGCG, CGGCGGGGGGGA, CGGTGGGGGGGG, TGGCGGGGGGGG]
		if (key.equals("CGGCGGGGGGGG")) {
			String [] expected = {"AGGCGGGGGGGG", "CGGCGGGGGGCG", "CGGCGGGGGGGA", "CGGTGGGGGGGG", "TGGCGGGGGGGG"};
			checkResult(key, values, expected);
		}
		// GGGGGGGGGGCC=[AGGGGGGGGGCC, CGGGGGGGGGCC, GGGGCGGGGGCC, GGGGGGGCGGCC, GGGGGGGGGCCC, GGGGGGGGGGAC, GGGGGGGGGGCA, GGGGGGGGGGCT, GGGGGGGGGTCC]
		if (key.equals("GGGGGGGGGGCC")) {
			String [] expected = {"AGGGGGGGGGCC", "CGGGGGGGGGCC", "GGGGCGGGGGCC", "GGGGGGGCGGCC", "GGGGGGGGGCCC", "GGGGGGGGGGAC", "GGGGGGGGGGCA", "GGGGGGGGGGCT", "GGGGGGGGGTCC"};
			checkResult(key, values, expected);
		}
		// GCGGCGGGGGGG=[ACGGCGGGGGGG, GAGGCGGGGGGG, GCGGCGGGGGCG, GCGGCGGGGGGA, GCGGCGGGGGGC, GTGGCGGGGGGG]
		if (key.equals("GCGGCGGGGGGG")) {
			String [] expected = {"ACGGCGGGGGGG", "GAGGCGGGGGGG", "GCGGCGGGGGCG", "GCGGCGGGGGGA", "GCGGCGGGGGGC", "GTGGCGGGGGGG"};
			checkResult(key, values, expected);
		}
		// AGGGGCGGGGGG=[AGAGGCGGGGGG, AGCGGCGGGGGG, AGGGGAGGGGGG, AGGGGCGGGGGA, AGGGGCGGTGGG, AGGGGTGGGGGG, AGTGGCGGGGGG, CGGGGCGGGGGG, TGGGGCGGGGGG]
		if (key.equals("AGGGGCGGGGGG")) {
			String [] expected = {"AGAGGCGGGGGG", "AGCGGCGGGGGG", "AGGGGAGGGGGG", "AGGGGCGGGGGA", "AGGGGCGGTGGG", "AGGGGTGGGGGG", "AGTGGCGGGGGG", "CGGGGCGGGGGG", "TGGGGCGGGGGG"};
			checkResult(key, values, expected);
		}
		// CCCCCCCCCCCC=[]
		if (key.equals("CCCCCCCCCCCC")) {
			Assert.assertTrue(values.size()==0);
		}
		// GGGGGGGGGGGG=[AGGGGGGGGGGG, CGGGGGGGGGGG, GAGGGGGGGGGG, GCGGGGGGGGGG, GGAGGGGGGGGG, GGCGGGGGGGGG, GGGAGGGGGGGG, GGGCGGGGGGGG, GGGGAGGGGGGG, GGGGCGGGGGGG, GGGGGAGGGGGG, GGGGGCGGGGGG, GGGGGGAGGGGG, GGGGGGCGGGGG, GGGGGGGAGGGG, GGGGGGGCGGGG, GGGGGGGGAGGG, GGGGGGGGCGGG, GGGGGGGGGAGG, GGGGGGGGGCGG, GGGGGGGGGGAG, GGGGGGGGGGCG, GGGGGGGGGGGA, GGGGGGGGGGGC, GGGGGGGGGGGT, GGGGGGGGGGTG, GGGGGGGGGTGG, GGGGGGGGTGGG, GGGGGGGTGGGG, GGGGGGTGGGGG, GGGGGTGGGGGG, GGGGTGGGGGGG, GGGTGGGGGGGG, GGTGGGGGGGGG, GTGGGGGGGGGG, TGGGGGGGGGGG]
		if (key.equals("GGGGGGGGGGGG")) {
			String [] expected = {"AGGGGGGGGGGG", "CGGGGGGGGGGG", "GAGGGGGGGGGG","GCGGGGGGGGGG", "GGAGGGGGGGGG", "GGCGGGGGGGGG", "GGGAGGGGGGGG", "GGGCGGGGGGGG", "GGGGAGGGGGGG", "GGGGCGGGGGGG", "GGGGGAGGGGGG", "GGGGGCGGGGGG", "GGGGGGAGGGGG", "GGGGGGCGGGGG", "GGGGGGGAGGGG", "GGGGGGGCGGGG", "GGGGGGGGAGGG", "GGGGGGGGCGGG", "GGGGGGGGGAGG", "GGGGGGGGGCGG", "GGGGGGGGGGAG", "GGGGGGGGGGCG", "GGGGGGGGGGGA", "GGGGGGGGGGGC", "GGGGGGGGGGGT", "GGGGGGGGGGTG", "GGGGGGGGGTGG", "GGGGGGGGTGGG", "GGGGGGGTGGGG", "GGGGGGTGGGGG", "GGGGGTGGGGGG", "GGGGTGGGGGGG", "GGGTGGGGGGGG", "GGTGGGGGGGGG", "GTGGGGGGGGGG", "TGGGGGGGGGGG"};
			checkResult(key, values, expected);
		}
		// GGTTTTTTTTTT=[CGTTTTTTTTTT, GCTTTTTTTTTT, GGATTTTTTTTT, GGGTTTTTTTTT]
		if (key.equals("GGTTTTTTTTTT")) {
			String [] expected = {"CGTTTTTTTTTT", "GCTTTTTTTTTT", "GGATTTTTTTTT", "GGGTTTTTTTTT"};
			checkResult(key, values, expected);
		}
		// CCGGGGGGGGGG=[ACGGGGGGGGGG, CAGGGGGGGGGG, CCCGGGGGGGGG, CCGGGGGGGGGA, CCGGGGGGGGGC, CTGGGGGGGGGG]
		if (key.equals("CCGGGGGGGGGG")) {
			String [] expected = {"ACGGGGGGGGGG", "CAGGGGGGGGGG", "CCCGGGGGGGGG", "CCGGGGGGGGGA", "CCGGGGGGGGGC", "CTGGGGGGGGGG"};
			checkResult(key, values, expected);
		}
		// AGGGGGGGGGGA=[AGGCGGGGGGGA, AGGGGAGGGGGA, AGGGGGCGGGGA, AGGGGGGCGGGA, AGGGGGGGAGGA, AGGGGGGGCGGA, AGGGGGGGGAGA, AGGGGGGGGCGA, AGGGGGGGGGAA, AGGGGGGGGGCA, AGGGGGGGGGGC, AGGGGGGGGGGT, CGGGGGGGGGGA, TGGGGGGGGGGA]
		if (key.equals("AGGGGGGGGGGA")) {
			String [] expected = {"AGGCGGGGGGGA", "AGGGGAGGGGGA", "AGGGGGCGGGGA", "AGGGGGGCGGGA", "AGGGGGGGAGGA", "AGGGGGGGCGGA", "AGGGGGGGGAGA", "AGGGGGGGGCGA", "AGGGGGGGGGAA", "AGGGGGGGGGCA", "AGGGGGGGGGGC", "AGGGGGGGGGGT", "CGGGGGGGGGGA", "TGGGGGGGGGGA"};
			checkResult(key, values, expected);
		}
		
	}
	
	private void checkResult (String key, List<String> values, String [] valuesExpected) {
		List<String> exp = new ArrayList<String>(Arrays.asList(valuesExpected));
		Assert.assertEquals(exp, values);
		
	}
	
	@Test(enabled=true)
	public void collapseBarcodesSmall() {
		ObjectCounter <String> barcodes = new ObjectCounter<String>();
		// primary barcode TEST1.
		barcodes.incrementByCount("TEST1", 10);
		barcodes.incrementByCount("TEST2", 9);
		barcodes.incrementByCount("TEST3", 8);
		// FEST3 primary barcode as TEST3 goes into TEST1.
		barcodes.incrementByCount("FEST3", 7);
		barcodes.incrementByCount("FEST2", 6);
		barcodes.incrementByCount("FEST4", 5);
		barcodes.incrementByCount("MEST3", 4);
		// MEST4 goes into FEST3, so MEST3 is a primary.
		barcodes.incrementByCount("MEST4", 3);
		barcodes.incrementByCount("MEST2", 2);
		// MEST1 goes into TEST1.
		barcodes.incrementByCount("MEST1", 1);
		
		MapBarcodesByEditDistance mapper = new  MapBarcodesByEditDistance(false, 4, 0);
		List<String>coreBarcodes  = barcodes.getKeysOrderedByCount(true);
		
		Map<String, List<String>> result = mapper.collapseBarcodes(coreBarcodes, barcodes, false, 1);
		Assert.assertEquals(3, result.size());
		
		Collection<String> r1 = result.get("TEST1");
		Assert.assertTrue(r1.contains("TEST2"));
		Assert.assertTrue(r1.contains("TEST3"));
		Assert.assertTrue(r1.contains("MEST1"));
		Assert.assertEquals(r1.size(),3);
		
		Collection<String> r2 = result.get("FEST3");
		Assert.assertTrue(r2.contains("FEST2"));
		Assert.assertTrue(r2.contains("FEST4"));
		Assert.assertTrue(r2.contains("MEST3"));
		Assert.assertEquals(r2.size(),3);
		
		
		Collection<String> r3 = result.get("MEST4");
		Assert.assertTrue(r3.contains("MEST2"));		
		Assert.assertEquals(r3.size(),1);
		
		Assert.assertNotNull(result);
		
	}
	
	 
	
}
