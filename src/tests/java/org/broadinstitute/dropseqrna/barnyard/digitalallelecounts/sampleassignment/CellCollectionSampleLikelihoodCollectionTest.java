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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IntervalList;
import junit.framework.Assert;
import org.apache.commons.lang.RandomStringUtils;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPUMIBasePileup;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPUMIBasePileupIterator;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SortOrder;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.CellCollectionSampleLikelihoodCollection;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.CellSampleLikelihoodCollection;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.SampleGenotypeProbabilities;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.SampleGenotypeProbabilitiesIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import org.testng.annotations.Test;
import picard.annotation.LocusFunction;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class CellCollectionSampleLikelihoodCollectionTest {

	private final static String rootDir="testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/sampleassignment/";
	
	private final File INPUT_BAM = new File(
			"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/smallTest_snpUMIPileUp_retagged.sam");
	private final File snpIntervalsFile = new File(
			"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/hek_cells_2snps.intervals");

	private final File likelihoodFile = new File(
			"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/sampleassignment/donor.assignments.txt");

	private final String cellBarcodeTag = "XC";
	private final String molBCTag = "XM";
	private final String snpTag = "YS";
	private final int readMQ = 10;

	private String GENE_NAME_TAG="gn";
	private String GENE_STRAND_TAG="gs";
	private String GENE_FUNCTION_TAG="gf";
	private StrandStrategy STRAND_STRATEGY = StrandStrategy.SENSE;
	private List<LocusFunction> LOCUS_FUNCTION_LIST=new ArrayList<LocusFunction>(Arrays.asList(LocusFunction.CODING, LocusFunction.UTR));

	@Test
	public void testParsing () {
		CellCollectionSampleLikelihoodCollection result = CellCollectionSampleLikelihoodCollection.parseFile(this.likelihoodFile);
		CellSampleLikelihoodCollection r = result.getLikelihoodCollection("CACGCAGGAGTG");
		Assert.assertEquals(9827, r.getNumSNPs());
		Assert.assertEquals(25077, r.getNumUMIs());
		// if you can infer these correctly, you probably did it right!
		Assert.assertEquals(1.0676E3, r.getBestSampleAssignment().getLogLikelihoodRatio(), 0.0001);
		Assert.assertEquals("WA13_P33_140530", r.getBestSampleAssignment().getSample());

		CellSampleLikelihoodCollection r2 = result.getLikelihoodCollection("GCAAGTCCAGCC");
		Assert.assertEquals(396, r2.getNumSNPs());
		Assert.assertEquals(453, r2.getNumUMIs());
		// if you can infer these correctly, you probably did it right!
		Assert.assertEquals(3.1779E1, r2.getBestSampleAssignment().getLogLikelihoodRatio(), 0.0001);
		Assert.assertEquals("Genea43_P12_150104", r2.getBestSampleAssignment().getSample());

	}
	
	@Test
	// a=read.table("testdata/org/broadinstitute/dropseq/private/barnyard/digitalallelecounts/sampleassignment/donor.assignments.txt", header=T, stringsAsFactors=F)
	// sort(a[a$cell=="CACGCAGGAGTG",7:dim(a)[2]], decreasing=T)
	public void testGetDonorsRankedByAssignmentLikelihood () {
		CellCollectionSampleLikelihoodCollection result = CellCollectionSampleLikelihoodCollection.parseFile(this.likelihoodFile);
		CellSampleLikelihoodCollection c = result.getLikelihoodCollection("CACGCAGGAGTG");		
		List<String> donorsByRank= c.getDonorsRankedByAssignmentLikelihood();
		// first 4 WA13_P33_140530 WA7_P33_140529 WA14_P21_131206 CSES15_P26_140611
		// last 4 MShef13_P16_14_0715 CHB8_P26_140723 HUES45_P21_131206 HES3_NKX2.1GFP_P10_13_0903
		List<String> expectedTop4 = Arrays.asList("WA13_P33_140530", "WA7_P33_140529", "WA14_P21_131206", "CSES15_P26_140611");
		List<String> expectedBottom4 = Arrays.asList("MShef13_P16_14_0715", "CHB8_P26_140723", "HUES45_P21_131206", "HES3_NKX2.1GFP_P10_13_0903");
		List<String> top4=donorsByRank.subList(0, 4);
		List<String> bot4 = donorsByRank.subList(125, 129);
		
		Assert.assertNotSame(expectedTop4, top4);
		Assert.assertNotSame(expectedBottom4, bot4);
		
	}

	@Test
	public void updateSNPUMICountsTest() {
		List<String> barcodes = new ArrayList<>();
		barcodes.add("ATCAGGGACAGA");
		barcodes.add("TTGCCTTACGCG");
		int editDistance=1;

		SampleGenotypeProbabilitiesIterator iter = prepIterFromFile(barcodes, editDistance);

		CellCollectionSampleLikelihoodCollection likelihoodCollection = new CellCollectionSampleLikelihoodCollection(Collections.EMPTY_LIST);
		while (iter.hasNext()) {
			SampleGenotypeProbabilities p = iter.next();
			likelihoodCollection.updateNumObservations(p);
		}
		// snp [HUMAN_1:76227022-76227022	+	HUMAN_1:76227022] cell [ATCAGGGACAGA] [G, A, G, A] [37, 8, 32, 8]
		// snp [HUMAN_1:76227022-76227022	+	HUMAN_1:76227022] cell [TTGCCTTACGCG] [G, A, A] [34, 37, 18]

		CellSampleLikelihoodCollection cslc = likelihoodCollection.getLikelihoodCollection("ATCAGGGACAGA");
		Assert.assertEquals(1, cslc.getNumSNPs());
		Assert.assertEquals(4, cslc.getNumUMIs());

		cslc = likelihoodCollection.getLikelihoodCollection("TTGCCTTACGCG");
		Assert.assertEquals(1, cslc.getNumSNPs());
		Assert.assertEquals(3, cslc.getNumUMIs());

	}

	@Test(enabled=false)
	public void testMemoryUsage () {
		List<String> barcodes = new ArrayList<>();
		barcodes.add("ATCAGGGACAGA");
		int editDistance=1;

		SampleGenotypeProbabilitiesIterator iter = prepIterFromFile(barcodes, editDistance);
		SampleGenotypeProbabilities p = iter.next();

		int numDonors = 150;
		int numCellBarcodes = 200000;
		int numSNPs = 200;
		char [] alleles = {'A', 'G'};
		Random random = new Random();

		List<String> donors = IntStream.range(1, numDonors+1).mapToObj(x -> RandomStringUtils.random(16, false, true)).collect(Collectors.toList());
		List<String> cellBarcodes = IntStream.range(1, numCellBarcodes+1).mapToObj(x -> RandomStringUtils.random(12, true, false)).collect(Collectors.toList());

		CellCollectionSampleLikelihoodCollection likelihoodCollection = new CellCollectionSampleLikelihoodCollection(donors);

		for (int s=0; s<numSNPs; s++) {
			if (s%1==0) System.out.println("SNP ["+ s +"|" + numSNPs +"]");
			for (int i=0; i<cellBarcodes.size(); i++) {
				String cellBarcode = cellBarcodes.get(i);
				SampleGenotypeProbabilities p2 = getPileUpWithNewCellBarcodeName(p, cellBarcode);
				for (int j=0; j<donors.size(); j++) {
					String sampleName = donors.get(j);
					char a1 = alleles[random.nextInt(alleles.length)];
					char a2 = alleles[random.nextInt(alleles.length)];
					likelihoodCollection.updatelikelihoods(p2, Collections.singletonList(sampleName), a1, a2, null, null);
				}
			} 
			System.gc();
		}
		Map<String, Double> result = likelihoodCollection.getFDRCorrectedPvalues();
		Assert.assertNotNull(result);

		// System.gc();
	    Runtime rt = Runtime.getRuntime();
	    long usedMB = (rt.totalMemory() - rt.freeMemory()) / 1024 / 1024;
	    System.out.println("Num SNPs [" +numSNPs +"] num donors [" + numDonors +"] num cells[" + numCellBarcodes+"] memory usage " + usedMB + "Mb");


	}

	private SampleGenotypeProbabilities getPileUpWithNewCellBarcodeName (final SampleGenotypeProbabilities p, final String cellBarcode) {
		SampleGenotypeProbabilities result = new SampleGenotypeProbabilities(p.getSNPInterval(), cellBarcode);
		for (SNPUMIBasePileup pile: p.getBackingPileups()) {
			SNPUMIBasePileup pile2 = new SNPUMIBasePileup(pile.getSNPInterval(), pile.getGene(), result.getCell(), pile.getMolecularBarcode());
			pile2.setBasesAndQualities(pile.getBases(), pile.getQualities());
			result.add(pile2);
		}
		result.toString();



		return result;
	}

	private SampleGenotypeProbabilitiesIterator prepIterFromFile (final List<String> barcodes, final int editDistance) {
		IntervalList snpIntervals = IntervalList.fromFile(snpIntervalsFile);
		
		
		
		SNPUMIBasePileupIterator sbpi = new SNPUMIBasePileupIterator(
				new SamHeaderAndIterator(this.INPUT_BAM), snpIntervals, GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG,
				LOCUS_FUNCTION_LIST, STRAND_STRATEGY, cellBarcodeTag,
				molBCTag, snpTag, "XF", readMQ, true, barcodes, null, SortOrder.SNP_CELL);

		boolean flag = sbpi.hasNext();
		final SAMSequenceDictionary dict = snpIntervals.getHeader().getSequenceDictionary();
		SampleGenotypeProbabilitiesIterator result = new SampleGenotypeProbabilitiesIterator(sbpi, dict, editDistance);
		return result;
	}

	
	/*  Disabled test as functionality is unused and disabled.
	@Test
	public void testAggregationStyle() {
		List<String> barcodes = new ArrayList<>();
		barcodes.add("ATCAGGGACAGA");
		barcodes.add("TTGCCTTACGCG");
		int editDistance=1;

		SampleGenotypeProbabilitiesIterator iter = prepIterFromFile(barcodes, editDistance);
		List<String> samples = Collections.EMPTY_LIST;
		
		CellCollectionSampleLikelihoodCollection likelihoodCollection = new CellCollectionSampleLikelihoodCollection(samples);

		while (iter.hasNext()) {
			CellCollectionSampleLikelihoodCollection current = new CellCollectionSampleLikelihoodCollection(samples);
			SampleGenotypeProbabilities p = iter.next();
			current.updateNumObservations(p);
			likelihoodCollection.add(current);
		}
		// snp [HUMAN_1:76227022-76227022	+	HUMAN_1:76227022] cell [ATCAGGGACAGA] [G, A, G, A] [37, 8, 32, 8]
		// snp [HUMAN_1:76227022-76227022	+	HUMAN_1:76227022] cell [TTGCCTTACGCG] [G, A, A] [34, 37, 18]

		CellSampleLikelihoodCollection cslc = likelihoodCollection.getLikelihoodCollection("ATCAGGGACAGA");
		Assert.assertEquals(1, cslc.getNumSNPs());
		Assert.assertEquals(4, cslc.getNumUMIs());

		cslc = likelihoodCollection.getLikelihoodCollection("TTGCCTTACGCG");
		Assert.assertEquals(1, cslc.getNumSNPs());
		Assert.assertEquals(3, cslc.getNumUMIs());


	}
	*/

}
