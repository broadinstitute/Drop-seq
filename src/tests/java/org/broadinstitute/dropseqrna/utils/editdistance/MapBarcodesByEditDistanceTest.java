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
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang3.RandomStringUtils;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.testng.annotations.Test;

import junit.framework.Assert;
import picard.util.TabbedInputParser;

public class MapBarcodesByEditDistanceTest {

	private static File testData = new File ("testdata/org/broadinstitute/transcriptome/utils/editdistance/inEditDistSmall.txt");

	private static File indelAnswerKey = new File ("testdata/org/broadinstitute/transcriptome/utils/editdistance/indel_barcode_repair_answer_key.txt");
	private static File repairedBC = new File ("testdata/org/broadinstitute/transcriptome/utils/editdistance/repairedBC.txt");
	private static File intendedBC = new File ("testdata/org/broadinstitute/transcriptome/utils/editdistance/potential_intendedBC.txt");

	@Test
	public void testFindIntendedIndelSequences() {
		// A set of data collapsed by R.
		TabbedInputParser parser = new TabbedInputParser(false, indelAnswerKey);

		//key: repaired barcode, value: intended barcode.
		Map<String, String> expectedResult = new HashMap<>();

		// skip the header
		String [] header = parser.next();
		while (parser.hasNext()) {
			String [] line = parser.next();
			expectedResult.put(line[0], line[1]);
		}


		// some abmiguous intended sequences that need to be in the starting data so the results work out correctly.
		List<String> intended=readFile(intendedBC);
		List<String> repaired = readFile(repairedBC);

		// now run the test!
		MapBarcodesByEditDistance mbed = new MapBarcodesByEditDistance(true);
		Map<String,String> result= mbed.findIntendedIndelSequences(repaired, intended, 1);

		// now test
		for (String repairedBC: result.keySet()) {
			String expectedIntended=expectedResult.get(repairedBC);
			if (expectedIntended==null)
				System.out.println("");
			String actualIntended =result.get(repairedBC);
			if (actualIntended==null)
				System.out.println("");
			Assert.assertNotNull(expectedIntended);
			Assert.assertNotNull(actualIntended);
			Assert.assertEquals(expectedIntended, actualIntended);
		}

		// check from the other direction
		for (String repairedBC: expectedResult.keySet()) {
			String expectedIntended=expectedResult.get(repairedBC);
			if (expectedIntended==null)
				System.out.println("");
			String actualIntended =result.get(repairedBC);
			if (actualIntended==null)
				System.out.println("");
			Assert.assertNotNull(expectedIntended);
			Assert.assertNotNull(actualIntended);
			Assert.assertEquals(expectedIntended, actualIntended);
		}

	}

	private List<String> readFile (final File f) {
		List<String> result = new ArrayList<>();
		TabbedInputParser parser = new TabbedInputParser(false, f);
		while (parser.hasNext()) {
			String [] line = parser.next();
			result.add(line[0]);
		}
		return result;
	}

	@Test
	public void testBottomUpCollapse() {
		MapBarcodesByEditDistance mbed=new MapBarcodesByEditDistance(true);
		ObjectCounter<String> barcodes = new ObjectCounter<>();
		// initialize with test data
		// unambiguous pair
		barcodes.incrementByCount("GTACAAAATATC", 1003);
		barcodes.incrementByCount("GTACAAAATATA", 4287);

		// unambiguous pair
		barcodes.incrementByCount("CCGCCGTTCGAA", 1016);
		barcodes.incrementByCount("CCGCAGTTCGAA", 3108);

		//ambiguous
		barcodes.incrementByCount("TAGAATCCCAAG", 20);
		barcodes.incrementByCount("TAGAATCACAAG", 39);
		barcodes.incrementByCount("TAGAATCGCAAG", 5199);

		//ambiguous
		barcodes.incrementByCount("CCGGAGACTATA", 20);
		barcodes.incrementByCount("CCGGAGACGATA", 23);
		barcodes.incrementByCount("CCGGAGACCATA", 2202);
		barcodes.incrementByCount("CCGGAGACAATA", 4165);

		Set<String> expectedAmbiguous = new HashSet<>(Arrays.asList("TAGAATCCCAAG", "CCGGAGACTATA", "CCGGAGACGATA"));


		// no neighbors
		barcodes.incrementByCount("ACTGTAGAAGGG", 172);

		// run and validate
		BottomUpCollapseResult result= mbed.bottomUpCollapse(barcodes, 1);

		// validate unambiguous
		String larger = result.getLargerRelatedBarcode("GTACAAAATATC");
		Assert.assertEquals("GTACAAAATATA", larger);
		larger = result.getLargerRelatedBarcode("CCGCCGTTCGAA");
		Assert.assertEquals("CCGCAGTTCGAA", larger);
		larger = result.getLargerRelatedBarcode("CCGGAGACCATA");
		Assert.assertEquals("CCGGAGACAATA", larger);
		larger = result.getLargerRelatedBarcode("TAGAATCACAAG");
		Assert.assertEquals("TAGAATCGCAAG", larger);

		// validate ambiguous
		Collection<String> ambiguous = result.getAmbiguousBarcodes();
		Assert.assertTrue(ambiguous.containsAll(expectedAmbiguous));
		Assert.assertTrue(expectedAmbiguous.containsAll(ambiguous));

		// validate no neighbor
		larger = result.getLargerRelatedBarcode("ACTGTAGAAGGG");
		Assert.assertNull(larger);

	}

	@Test(enabled=false)
	public void testBottomUpSpeed () {
		MapBarcodesByEditDistance mbed=new MapBarcodesByEditDistance(true, 1, 10000);
		ObjectCounter<String> barcodes = getRandomBarcodes(12, 500000);
		BottomUpCollapseResult result= mbed.bottomUpCollapse(barcodes, 1);
		Assert.assertNotNull(result);
	}

	private ObjectCounter<String> getRandomBarcodes (final int barcodeLength, final int numBarcodes) {
		List<String> barcodes = getRandomBarcodesAsList(barcodeLength, numBarcodes);
		ObjectCounter<String> b = new ObjectCounter<>();
		barcodes.stream().forEach(x -> b.increment(x));
		return b;
	}

	public static List<String> getRandomBarcodesAsList (final int numBases, final int numBarcodes) {
		char [] bases = {'A', 'C', 'G', 'T', 'N'};
		List<String> result = new ArrayList<>(numBarcodes);
		for (int i=0; i<numBarcodes; i++)
			result.add(RandomStringUtils.random(numBases, bases));
		return (result);
	}

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



	private void checkResult (final String key, final List<String> values) {
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
		if (key.equals("CCCCCCCCCCCC"))
			Assert.assertTrue(values.size()==0);
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

	private void checkResult (final String key, final List<String> values, final String [] valuesExpected) {
		List<String> exp = new ArrayList<>(Arrays.asList(valuesExpected));
		Assert.assertEquals(exp, values);

	}

	@Test(enabled=true)
	public void collapseBarcodesSmall() {
		ObjectCounter <String> barcodes = new ObjectCounter<>();
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

	@Test
	public void testCollapseBarcodes () {
		MapBarcodesByEditDistance m = new MapBarcodesByEditDistance(false);
		ObjectCounter<String> barcodes= new ObjectCounter<>();

		barcodes.incrementByCount("AAACCCGGGTTT", 20);  // winner
		barcodes.incrementByCount("AAACCCGGAGTT", 10);  // indel 1.
		barcodes.incrementByCount("AAACCCGGGTAT", 8);  // sub 1
		barcodes.incrementByCount("AAACCCGGGCTT", 6);  // sub 1.
		// group 2.
		barcodes.incrementByCount("AAACGGGAGGTA", 2);
		barcodes.incrementByCount("GTAGACTAGGTG", 2);
		barcodes.incrementByCount("TGCAGGTGGCCG",2);
		barcodes.incrementByCount("CCGTGCGTCACA", 2);
		barcodes.incrementByCount("GGGCCCTATCAT", 2);
		barcodes.incrementByCount("CACTGCCGTTGG", 2);

		Map<String, List<String>> result = m.collapseBarcodes(barcodes, true, 1);
		List<String> expected = Arrays.asList("AAACCCGGAGTT", "AAACCCGGGCTT", "AAACCCGGGTAT");
		List<String> actual = result.get("AAACCCGGGTTT");
		Assert.assertEquals(expected, actual);

		Assert.assertNull(result.get("AAACCCGGAGTT"));
		Assert.assertNull(result.get("AAACCCGGGCTT"));
		Assert.assertNull(result.get("AAACCCGGGTAT"));



	}

	@Test
	public void testCollapseAndMergeBarcodes () {
		MapBarcodesByEditDistance m = new MapBarcodesByEditDistance(false);
		ObjectCounter<String> barcodes= new ObjectCounter<>();

		barcodes.incrementByCount("AAACCCGGGTTT", 20);  // winner
		barcodes.incrementByCount("AAACCCGGAGTT", 10);  // indel 1.
		barcodes.incrementByCount("AAACCCGGGTAT", 8);  // sub 1
		barcodes.incrementByCount("AAACCCGGGCTT", 6);  // sub 1.
		// group 2.
		barcodes.incrementByCount("AAACGGGAGGTA", 2);
		barcodes.incrementByCount("GTAGACTAGGTG", 2);
		barcodes.incrementByCount("TGCAGGTGGCCG",2);
		barcodes.incrementByCount("CCGTGCGTCACA", 2);
		barcodes.incrementByCount("GGGCCCTATCAT", 2);
		barcodes.incrementByCount("CACTGCCGTTGG", 2);

		ObjectCounter<String> result = m.collapseAndMergeBarcodes(barcodes, true, 1);
		Assert.assertEquals(44, result.getCountForKey("AAACCCGGGTTT"));
		Assert.assertEquals(0, result.getCountForKey("AAACCCGGAGTT"));
		Assert.assertEquals(0, result.getCountForKey("AAACCCGGGTAT"));
		Assert.assertEquals(0, result.getCountForKey("AAACCCGGGCTT"));

		Assert.assertEquals(2, result.getCountForKey("AAACGGGAGGTA"));
		Assert.assertEquals(2, result.getCountForKey("GTAGACTAGGTG"));
		Assert.assertEquals(2, result.getCountForKey("TGCAGGTGGCCG"));

	}

	@Test
	public void testFindEditDistanceThreshold () {
		// Group 1:
		// z=c("AAACCCGGGTTT", "AAACCCGGGTTA", "AAACCCGGGTAT", "AAACCCGGGCTT")

		// Group 2:
		// z2=c("AAACGGGAGGTA","GTAGACTAGGTG","TGCAGGTGGCCG","CCGTGCGTCACA","GGGCCCTATCAT","CACTGCCGTTGG")

		MapBarcodesByEditDistance m = new MapBarcodesByEditDistance(true);
		Set<String> barcodes= new HashSet<>();

		barcodes.add("AAACCCGGAGTT");  // indel 1.
		barcodes.add("AAACCCGGGTAT");  // sub 1
		barcodes.add("AAACCCGGGCTT");  // sub 1.
		// group 2.
		barcodes.add("AAACGGGAGGTA");
		barcodes.add("GTAGACTAGGTG");
		barcodes.add("TGCAGGTGGCCG");
		barcodes.add("CCGTGCGTCACA");
		barcodes.add("GGGCCCTATCAT");
		barcodes.add("CACTGCCGTTGG");


		int threshold = m.findEditDistanceThreshold("AAACCCGGGTTT", barcodes, true, 1, 3);
		Assert.assertEquals(1, threshold);

		threshold = m.findEditDistanceThreshold("AAACCCGGGTTT", barcodes, false, 1, 3);
		Assert.assertEquals(2, threshold);

	}
	@Test (enabled=true)
	public void testCollapseBarcodesAdaptive () {
		// Group 1:
		// z=c("AAACCCGGGTTT", "AAACCCGGGTTA", "AAACCCGGGTAT", "AAACCCGGGCTT")

		// Group 2:
		// z2=c("AAACGGGAGGTA","GTAGACTAGGTG","TGCAGGTGGCCG","CCGTGCGTCACA","GGGCCCTATCAT","CACTGCCGTTGG")

		MapBarcodesByEditDistance m = new MapBarcodesByEditDistance(true,2,5);
		ObjectCounter<String> barcodes= new ObjectCounter<>();

		barcodes.incrementByCount("AAACCCGGGTTT", 20);  // the winner.
		barcodes.incrementByCount("AAACCCGGAGTT", 18);  // indel 1.
		barcodes.incrementByCount("AAACCCGGGTAT", 15);  // sub 1
		barcodes.incrementByCount("AAACCCGGGCTT", 13);  // sub 1.

		barcodes.incrementByCount("AAACGGGAGGTA", 2);   //
		barcodes.incrementByCount("GTAGACTAGGTG", 2);
		barcodes.incrementByCount("TGCAGGTGGCCG", 2);
		barcodes.incrementByCount("CCGTGCGTCACA", 2);
		barcodes.incrementByCount("GGGCCCTATCAT", 2);
		barcodes.incrementByCount("CACTGCCGTTGG", 2);


		MapBarcodesByEditDistance.AdaptiveMappingResult r = m.collapseBarcodesAdaptive(barcodes, true, 3, 1, 3);
		Map<String, List<String>> collapse = r.getBarcodeCollapseResult();

		Assert.assertTrue(collapse!=null);
		// we expect the winner to own the next 3 barcodes.
		// sorted alphabetically.
		List<String> expected = Arrays.asList("AAACCCGGAGTT", "AAACCCGGGCTT", "AAACCCGGGTAT");
		List<String> actual = collapse.get("AAACCCGGGTTT");
		Assert.assertEquals(expected, actual);

		// these were merged.
		Assert.assertNull(collapse.get("AAACCCGGAGTT"));
		Assert.assertNull(collapse.get("AAACCCGGGCTT"));
		Assert.assertNull(collapse.get("AAACCCGGGTAT"));

		// but those far away barcodes aren't merged.
		Assert.assertEquals(0, collapse.get("AAACGGGAGGTA").size());
		Assert.assertEquals(0, collapse.get("CACTGCCGTTGG").size());


	}



}




