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
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang3.RandomStringUtils;
import org.broadinstitute.dropseqrna.metrics.UmiSharingMetrics;
import org.broadinstitute.dropseqrna.metrics.umisharing.ParentEditDistanceMatcher;
import org.broadinstitute.dropseqrna.metrics.umisharing.ParentEditDistanceMatcher.TagValues;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;
import junit.framework.Assert;
import picard.util.TabbedInputParser;

public class MapBarcodesByEditDistanceTest {

	private static File testData = new File ("testdata/org/broadinstitute/transcriptome/utils/editdistance/inEditDistSmall.txt");

	private static File indelAnswerKey = new File ("testdata/org/broadinstitute/transcriptome/utils/editdistance/indel_barcode_repair_answer_key.txt");
	private static File repairedBC = new File ("testdata/org/broadinstitute/transcriptome/utils/editdistance/repairedBC.txt");
	private static File intendedBC = new File ("testdata/org/broadinstitute/transcriptome/utils/editdistance/potential_intendedBC.txt");

	
	private static File UmiSharingData = new File("testdata/org/broadinstitute/dropseq/metrics/compute_umi_sharing.unmapped.sam");
	private static File UmiSharingResultsED0 = new File ("testdata/org/broadinstitute/transcriptome/utils/editdistance/umi_test_data.merged_barcodes_ed0.txt");
	private static File UmiSharingResultsED1 = new File ("testdata/org/broadinstitute/transcriptome/utils/editdistance/umi_test_data.merged_barcodes_ed1.txt");
	
	// there's also testing of UMISharing in ComputeUMISharingTest...
	
	//broad/mccarroll/software/dropseq/jn_branch/BamTagOfTagCounts I=compute_umi_sharing.unmapped.sam PRIMARY_TAG=rb SECONDARY_TAG=XM FILTER_PCR_DUPLICATES=false O=compute_umi_sharing.unmapped.rb:XM.histogram.txt MINIMUM_MAPPING_QUALITY=0

	@Test
	public void testUMISharingCollapseED0() {
		testUMISharingCollapse(this.UmiSharingData, 0, this.UmiSharingResultsED0);
	}
	
	@Test
	public void testUMISharingCollapseED1() {
		testUMISharingCollapse(this.UmiSharingData, 1, this.UmiSharingResultsED1);
	}
	
	public void testUMISharingCollapse (File input, int editDistance, File resultFile) {		
		ParentEditDistanceMatcher matcher = new ParentEditDistanceMatcher(Collections.singletonList("XM"), Collections.singletonList(editDistance), false, 1);
		double sharingThreshold = 0.8d;
		UMISharingData data = readUMISharingData(UmiSharingData, "rb", matcher);
		MapBarcodesByEditDistance mbed = new MapBarcodesByEditDistance(false);		
		FindSimilarEntitiesResult<String, UmiSharingMetrics> result =mbed.collapseBarcodesByUmiSharing(data.barcodes, matcher, sharingThreshold, data.umisPerBarcode);
		Set<UmiSharingMetrics> metricsList = result.getCollapseMetric();
		// filter the metrics to those with at least some sharing, which is what the R code does.
		metricsList = metricsList.stream().filter(x -> x.FRAC_SHARED>=sharingThreshold).collect(Collectors.toSet());		
		FindSimilarEntitiesResult<String, UmiSharingMetrics> expectedResult = readExpectedResult(resultFile);
		
		// all the expected sharing is there.
		for (UmiSharingMetrics m : expectedResult.getCollapseMetric()) {						
			Assert.assertTrue(metricsList.contains(m));
		}
		
		Map<String, List<String>> mapping = result.getEntityMap();
		Map<String, List<String>> expectedMapping = expectedResult.getEntityMap();
		Assert.assertTrue(mapping.keySet().containsAll(expectedMapping.keySet()));
		
		for (String key: expectedMapping.keySet()) {
			Set<String> values = new HashSet<String> (mapping.get(key));
			Set<String> expectedValues = new HashSet<String> (expectedResult.getEntityMap().get(key));
			Assert.assertEquals(expectedValues, values);
		}		
		Assert.assertNotNull(result);		
	}
	
	private FindSimilarEntitiesResult<String, UmiSharingMetrics> readExpectedResult (File expected) {
		TabbedInputParser p = new TabbedInputParser(true, expected);
		String [] header = p.next();
		Set<UmiSharingMetrics> collapseMetric = new HashSet<UmiSharingMetrics>();
		
		Map<String, List<String>> mapping = new HashMap<String,List<String>>();
		FindSimilarEntitiesResult<String, UmiSharingMetrics> result = new FindSimilarEntitiesResult<>();
		
		while (p.hasNext()) {
			String [] line = p.next();
			String s= line[0];
			UmiSharingMetrics m = new UmiSharingMetrics();
			m.PARENT=line[0];
			m.CHILD=line[1];
			m.NUM_PARENT= Integer.parseInt(line[2]);
			m.NUM_CHILD=Integer.parseInt(line[3]);
			m.NUM_SHARED=Integer.parseInt(line[4]);
			m.FRAC_SHARED=Double.parseDouble(line[6]);
			result.addMetrics(m);
			result.addMapping(m.PARENT, m.CHILD);			
		}		
		result.addMapping(mapping);
		return result;
	}
	
	private UMISharingData readUMISharingData (File f, String cellBarcode, ParentEditDistanceMatcher matcher) {
		Map<String, Set<TagValues>> umisPerBarcode = new HashMap<String, Set<TagValues>>();
		ObjectCounter<String> barcodes = new ObjectCounter<>();
		
		SamReader reader = SamReaderFactory.makeDefault().open(UmiSharingData);
		SAMRecordIterator iter = reader.iterator();
		while (iter.hasNext()) {
			SAMRecord r = iter.next();
			String barcode = r.getStringAttribute(cellBarcode);
			TagValues v = matcher.getValues(r);
			Set<TagValues> s = umisPerBarcode.get(barcode);
			if (s==null) {
				s = new HashSet<TagValues>();
				umisPerBarcode.put(barcode, s);
			}			
			if (!s.contains(v))
				barcodes.increment(barcode);
			s.add(v);
		}
		UMISharingData result = new UMISharingData();
		result.barcodes=barcodes;
		result.umisPerBarcode=umisPerBarcode;
		return result;
	}
	
	private class UMISharingData {
		Map<String, Set<TagValues>> umisPerBarcode;
		ObjectCounter<String> barcodes;
	}

	@Test
	public void testFindRelatedBarcodesByMutationalCollapse1() {
		String coreBC = "ACGTG";
		// ACGTA->AAGTA=1 (accept)
		// ACGTA->AATTA=2 & AAGTA->AATTA=1 (accept)
		// ACGTA->AATTA=3 & AAGTA->AATTA=1 (accept)
		// ACGTA->TATTC=4 & AATTA->TATTC=2 (reject)
		List<String> relatedBarcodes=Arrays.asList("ACGTA", "AAGTA", "AATTA", "TATTC");
		MapBarcodesByEditDistance m = new MapBarcodesByEditDistance(true);
		FindSimilarEntitiesByMutationalCollapse func = new FindSimilarEntitiesByMutationalCollapse(m, false, 5, 1);
		List<String> actual = func.find(coreBC, relatedBarcodes, null).getEntityMap().get(coreBC);
		Collections.sort(actual);
		List<String> expected = Arrays.asList("ACGTA", "AAGTA", "AATTA");
		Collections.sort(expected);
		Assert.assertEquals(expected, actual);
	}
	
	@Test
	public void testFindRelatedBarcodesByMutationalCollapseED2() {
		// AAAAAA-> AAAGGA=1 (accept)
		// AAAAAA-> ATTGGA=4 & AAAGGA-> ATTGGA=2 (accept)
		// AAAAAA-> CTTGGG=6 & ATTGGA->CTTGGG=2 (accept)
		// AAAAAA-> CTTTTT=6 & CTTGGG->CTTTTT=3 (reject)
		String coreBC = "AAAAAA";
		List<String> relatedBarcodes=Arrays.asList("AAAGGA", "ATTGGA", "CTTGGG", "CTTTTT");
		
		MapBarcodesByEditDistance m = new MapBarcodesByEditDistance(true);
		FindSimilarEntitiesByMutationalCollapse func = new FindSimilarEntitiesByMutationalCollapse(m, false, 6, 1);
		Set<String> actual1 = new HashSet<String>(func.find(coreBC, relatedBarcodes, null).getEntityMap().get(coreBC));		
		Set<String> expected1 = new HashSet<String>(Arrays.asList());
		Assert.assertEquals(expected1, actual1);
		
		func = new FindSimilarEntitiesByMutationalCollapse(m, false, 6, 2);
		Set<String> actual = new HashSet<String>(func.find(coreBC, relatedBarcodes, null).getEntityMap().get(coreBC));			
		Set<String> expected = new HashSet<String>(Arrays.asList("AAAGGA", "ATTGGA", "CTTGGG"));
		Assert.assertEquals(expected, actual);				
		
	}


	@Test(enabled=true)
	public void testFindRelatedBarcodesByMutationalCollapse2() {
		MapBarcodesByEditDistance m = new MapBarcodesByEditDistance(true);

		// core one:
		// ACGTAAT->AAGTAAT=1 (accept)
		// ACGTAAT->AATTAAT=2 & AAGTAAT->AATTAAT=1 (accept)
		// ACGTAAT->AATTATT=3 & AAGTAAT->AATTAAT=1 (accept)
		// ACGTAAT->TATTATG=5 & AATTATT->TATTATG=2 (reject)
		List<String> relatedBarcodesOne=Arrays.asList("AAGTAAT", "AATTAAT", "AATTATT", "TATTATG");
		String coreBC="ACGTAAT";
		FindSimilarEntitiesByMutationalCollapse func = new FindSimilarEntitiesByMutationalCollapse(m, false, 5, 1);
		Set<String> actual = new HashSet<String> (func.find(coreBC, relatedBarcodesOne, null).getEntityMap().get(coreBC));						
		Set<String> expectedOne = new HashSet<String> (Arrays.asList("AAGTAAT", "AATTAAT", "AATTATT"));
		Assert.assertEquals(expectedOne, actual);

		// core two:
		// GGCCATG->GACCATG=1 (accept)
		// GGCCATG->GACCATA=2 & GACCATG->GACCATA=1 (accept)
		// GGCCATG->AACCATA=3 & GACCATA->AACCATA=1 (accept)
		// GGCCATG->AATCATA=4 & AACCATA->AATCATA=1 (accept)
		// GGCCATG->ACTCATA=4 & AGCCATA->ACTCATA=2 (reject)
		List<String> relatedBarcodesTwo=Arrays.asList("GACCATG", "GACCATA", "AACCATA", "AATCATA", "ACTCATA");
		coreBC="GGCCATG";		
		actual = new HashSet<String>(func.find(coreBC, relatedBarcodesTwo, null).getEntityMap().get(coreBC));		
		Set<String> expectedTwo = new HashSet<String> (Arrays.asList("GACCATG", "GACCATA", "AACCATA", "AATCATA"));
		Assert.assertEquals(expectedTwo, actual);

		// core three:
		// TTGAACA->TTGAACT=1 (accept)
		// TTGAACA->ATGAACT=2 & TTGAACT->ATGAACT=1 (accept)
		// TTGAACA->ATGTACT=3 & ATGAACT->ATGTACT=1 (accept)
		// TTGAACA->AGGTACT=4 & ATGTACT->AGGTACT=1 (accept)
		// TTGAACA->TGCTACT=4 & AGGTACT->TGCTACT=2 (reject)
		List<String> relatedBarcodesThree=Arrays.asList("TTGAACT", "ATGAACT", "ATGTACT", "AGGTACT", "TGCTACT");
		coreBC="TTGAACA";
		actual = new HashSet<String> (func.find(coreBC, relatedBarcodesThree, null).getEntityMap().get(coreBC));		
		Set<String> expectedThree = new HashSet<String> (Arrays.asList("TTGAACT", "ATGAACT", "ATGTACT", "AGGTACT"));
		Assert.assertEquals(expectedThree, actual);
		
		// test all of them.
		// plot (hclust(stringdistmatrix(c("ACGTAAT", "GGCCATG", "TTGAACA", "AAGTAAT", "AATTAAT", "AATTATT", "TATTATG", "GACCATG", "GACCATA", "AACCATA", "AATCATA", "ACTCATA", "TTGAACT", "ATGAACT", "ATGTACT", "AGGTACT", "TGCTACT"), method="hamming", useNames="strings")))

		List<String> relatedBarcodesAll = new ArrayList<String>();
		relatedBarcodesAll.addAll(relatedBarcodesOne);
		relatedBarcodesAll.addAll(relatedBarcodesTwo);
		relatedBarcodesAll.addAll(relatedBarcodesThree);
		
		coreBC="ACGTAAT";
		actual = new HashSet<String> (func.find(coreBC, relatedBarcodesOne, null).getEntityMap().get(coreBC));		
		Assert.assertEquals(expectedOne, actual);

		coreBC="GGCCATG";
		actual = new HashSet<String> (func.find(coreBC, relatedBarcodesTwo, null).getEntityMap().get(coreBC));		
		Assert.assertEquals(expectedTwo, actual);
		
		coreBC="TTGAACA";
		actual = new HashSet<String> (func.find(coreBC, relatedBarcodesThree, null).getEntityMap().get(coreBC));
		Assert.assertEquals(expectedThree, actual);
		
	}

	@Test 
	public void testCollapseBarcodesByMutationalCollapse () {
		
		MapBarcodesByEditDistance m = new MapBarcodesByEditDistance(true);

		// core one:
		// ACGTAAT->AAGTAAT=1 (accept)
		// ACGTAAT->AATTAAT=2 & AAGTAAT->AATTAAT=1 (accept)
		// ACGTAAT->AATTATT=3 & AAGTAAT->AATTAAT=1 (accept)
		// ACGTAAT->TATTATG=5 & AATTATT->TATTATG=2 (reject)
		
		// core two:
		// GGCCATG->GACCATG=1 (accept)
		// GGCCATG->GACCATA=2 & GACCATG->GACCATA=1 (accept)
		// GGCCATG->AACCATA=3 & GACCATA->AACCATA=1 (accept)
		// GGCCATG->AATCATA=4 & AACCATA->AATCATA=1 (accept)
		// GGCCATG->ACTCATA=4 & AGCCATA->ACTCATA=2 (reject)
		
		// core three:
		// TTGAACA->TTGAACT=1 (accept)
		// TTGAACA->ATGAACT=2 & TTGAACT->ATGAACT=1 (accept)
		// TTGAACA->ATGTACT=3 & ATGAACT->ATGTACT=1 (accept)
		// TTGAACA->AGGTACT=4 & ATGTACT->AGGTACT=1 (accept)
		// TTGAACA->TGCTACT=4 & AGGTACT->TGCTACT=2 (reject)
				
		// {TGCTACT=[], ACTCATA=[], GGCCATG=[AACCATA, AATCATA, GACCATA, GACCATG], TATTATG=[], ACGTAAT=[AAGTAAT, AATTAAT, AATTATT], TTGAACA=[AGGTACT, ATGAACT, ATGTACT, TTGAACT]}
		
		// the data.
		ObjectCounter<String> barcodes = new ObjectCounter<>();
		barcodes.incrementByCount("ACGTAAT", 90);
		
		barcodes.incrementByCount("AAGTAAT", 10);
		barcodes.incrementByCount("AATTAAT", 9);
		barcodes.incrementByCount("AATTATT", 7);
		barcodes.incrementByCount("TATTATG", 5);
		
		barcodes.incrementByCount("GGCCATG", 95);
		
		barcodes.incrementByCount("GACCATG", 10);
		barcodes.incrementByCount("GACCATA", 9);
		barcodes.incrementByCount("AACCATA", 7);
		barcodes.incrementByCount("AATCATA", 5);
		barcodes.incrementByCount("ACTCATA", 4);
				
		barcodes.incrementByCount("TTGAACA", 78);
		
		barcodes.incrementByCount("TTGAACT", 10);
		barcodes.incrementByCount("ATGAACT", 9);
		barcodes.incrementByCount("ATGTACT", 7);
		barcodes.incrementByCount("AGGTACT", 5);
		barcodes.incrementByCount("TGCTACT", 5);
		
		Map<String, List<String>> result = m.collapseBarcodesByMutationalCollapse (barcodes, false, 5, 3, 1);
		
		// expected results.
		String coreBCOne="ACGTAAT";
		String coreBCTwo="GGCCATG";
		String coreBCThree="TTGAACA";
		
		Set<String> expectedOne = new HashSet<>(Arrays.asList("AAGTAAT", "AATTAAT", "AATTATT"));
		Set<String> expectedTwo = new HashSet<>(Arrays.asList("GACCATG", "GACCATA", "AACCATA", "AATCATA"));
		Set<String> expectedThree = new HashSet<>(Arrays.asList("TTGAACT", "ATGAACT", "ATGTACT", "AGGTACT"));
		
		Assert.assertTrue(result.containsKey(coreBCOne));
		Assert.assertTrue(result.containsKey(coreBCTwo));
		Assert.assertTrue(result.containsKey(coreBCThree));
		
		// 3 cores + 3 unmerged.
		
		Assert.assertEquals(6, result.keySet().size());
		Assert.assertEquals(expectedOne, new HashSet<String>(result.get(coreBCOne)));
		Assert.assertEquals(expectedTwo, new HashSet<String>(result.get(coreBCTwo)));
		Assert.assertEquals(expectedThree, new HashSet<String>(result.get(coreBCThree)));
	}
	
	
	
	@Test(enabled=false)
	//TODO: re-enable this test when we agree on the corrrect answer.
	public void testCollapseBarcodesByMutationalCollapseLargeInput() {
		File dataFile = new File ("testdata/org/broadinstitute/transcriptome/utils/editdistance/mutational_collapse_testdata.txt.gz");
		// File dataFile = new File ("testdata/org/broadinstitute/transcriptome/utils/editdistance/170330_pSPBN_GFP_v9_v2_B19EnvA_15P_BCpooled_day5_Final_RVg_barcode.counts.txt.gz");
		File resultFile = new File ("testdata/org/broadinstitute/transcriptome/utils/editdistance/mutational_collapse_testdata.result.txt.gz");
		int minCount=3;
		
		TabbedInputParser dataParser = new TabbedInputParser(false, dataFile);				
		ObjectCounter<String> data = new ObjectCounter<>();		
		dataParser.next();
		while (dataParser.hasNext()) {
			String [] line = dataParser.next();
			data.incrementByCount(line[0], Integer.parseInt(line[1]));
		}
		dataParser.close();
		
		TabbedInputParser resultParser = new TabbedInputParser(false, resultFile);
		ObjectCounter<String> expectedRawCount = new ObjectCounter<>();	
		ObjectCounter<String> expectedAggregateCount = new ObjectCounter<>();
		resultParser.next();
		while (resultParser.hasNext()) {
			String [] line = resultParser.next();
			expectedRawCount.incrementByCount(line[0], Integer.parseInt(line[1]));
			expectedAggregateCount.incrementByCount(line[0], Integer.parseInt(line[5]));
		}
		resultParser.close();
		
		MapBarcodesByEditDistance m = new MapBarcodesByEditDistance(true, 4, 100);					
		Map<String, List<String>> result = m.collapseBarcodesByMutationalCollapse (data, false, 5, minCount, 1);
		
		ObjectCounter<String> r2 = aggregateCounts(data, result);
		writeOutput (data, r2, result, m);
						
		// Assert.assertEquals(r2.getKeys().size(), expected.getKeys().size());
		int numErrors=0;
		// compare two results.
		for (String key: expectedRawCount.getKeys()) {
			int e = expectedAggregateCount.getCountForKey(key);
			int a = r2.getCountForKey(key);
			
			if (e!=a) {				
//				List<String> children = result.get(key);
//				System.out.println("Mapping doesn't match for parent: " + key +" children "+ children.toString());
				numErrors++;
			}
			// Assert.assertEquals(e, a);
		}
		Assert.assertEquals(0, numErrors);
		Assert.assertNotNull(result);
						
	}
	
	private void writeOutput (ObjectCounter<String> data, ObjectCounter<String> aggregateCounts, Map<String, List<String>> mapping, MapBarcodesByEditDistance m) {				
		PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(new File ("/downloads/arpy/mutational_collapse_testdata.result_jim.txt")));			
		String [] header = {"sequence",  "counts", "parent",  "edist",  "fam_seqs", "fam_counts"};
		out.println(StringUtils.join(header, "\t"));
		
		for (String parentSeq: mapping.keySet()) {		
			if (parentSeq.equals("ACCTAGATTCAGTTAGCCAG"))
					System.out.println();
			List<String> sequences = mapping.get(parentSeq);														
			int famSeqs=sequences.size()+1;
			String [] line = {parentSeq, Integer.toString(data.getCountForKey(parentSeq)), parentSeq, "0", Integer.toString(famSeqs) ,Integer.toString(aggregateCounts.getCountForKey(parentSeq))};
			out.println(StringUtils.join(line, "\t"));
															
			for (String v: sequences) {
				int ed = HammingDistance.getHammingDistance(parentSeq, v);
				// for merged results, the family seqs size is always 1 and the fam counts is always 0.
				String [] line2 = {v, Integer.toString(data.getCountForKey(v)), parentSeq, Integer.toString(ed), "1", "0"};
				out.println(StringUtils.join(line2, "\t"));																				
			}			
		}
				
		out.close();
	}
	
	private ObjectCounter<String> aggregateCounts (ObjectCounter<String> data, Map<String, List<String>> mapping) {
		ObjectCounter<String> result = new ObjectCounter<>();
		// build out all the new counts for barcodes that have merged results.
		Set<String> mergedBarcodes=new HashSet<String>();
		
		for (String key: mapping.keySet()) {
			List<String> values = mapping.get(key);
			int totalCount=data.getCountForKey(key);			
			mergedBarcodes.addAll(values);
			totalCount += values.stream().mapToInt(x-> data.getCountForKey(x)).sum();						
			result.incrementByCount(key, totalCount);
			mergedBarcodes.add(key);			
		}
		
		// build out the singletons that were not merged by finding the non-merged data keys.
		for (String k: data.getKeys()) 
			if (!mergedBarcodes.contains(k)) {
				result.incrementByCount(k, data.getCountForKey(k));
			}
		
		
		return result;
		
	}

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

	@Test
	public void testBottomUpCollapseThreaded() {
		MapBarcodesByEditDistance mbed=new MapBarcodesByEditDistance(true, 100, 2);
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

		MapBarcodesByEditDistance mapper = new  MapBarcodesByEditDistance(false, 1, 0);
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
		FindSimilarEntitiesByAdaptiveEditDistance func = new FindSimilarEntitiesByAdaptiveEditDistance(m, true, 1, 1, 3);
		
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


		int threshold = func.findEditDistanceThreshold("AAACCCGGGTTT", barcodes);
		Assert.assertEquals(1, threshold);

		func = new FindSimilarEntitiesByAdaptiveEditDistance(m, false, 1, 1, 3);		
		threshold = func.findEditDistanceThreshold("AAACCCGGGTTT", barcodes);
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

	@Test(expectedExceptions=IllegalArgumentException.class)
	public void testBarcodesSameLength () {
		MapBarcodesByEditDistance m = new MapBarcodesByEditDistance(true,2,5);
		ObjectCounter <String> barcodes = new ObjectCounter<>();
		barcodes.increment("FOO");
		barcodes.increment("TOOLONG");
		m.bottomUpCollapse (barcodes, 1);
	}

}




