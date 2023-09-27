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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import picard.annotation.LocusFunction;

public class SNPUMIBasePileupIteratorTest {

	// See smallTest_snpUMIPileUp.summary.txt for a listing of the bases/qualities/cell/molecular barcodes for all 15 reads
	private final File smallBAMFile = new File(
			"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/smallTest_snpUMIPileUp.retagged.bam");
	private final File cellBCFile = new File(
			"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/hek_cells_cell_barcodes.txt");
	private final File snpIntervalsFile = new File(
			"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/hek_cells_2snps.intervals");
	private final String cellBarcodeTag = "XC";
	private final String molBCTag = "XM";
	private final String snpTag = "YS";
	private final String functionTag="XF";
	private final int readMQ = 10;
	private final boolean assignReadsToAllGenes = true;

	private String GENE_NAME_TAG="gn";
	private String GENE_STRAND_TAG="gs";
	private String GENE_FUNCTION_TAG="gf";
	private StrandStrategy STRAND_STRATEGY = StrandStrategy.SENSE;
	private List<LocusFunction> LOCUS_FUNCTION_LIST=new ArrayList<LocusFunction>(Arrays.asList(LocusFunction.CODING, LocusFunction.UTR));

	@Test
	public void testPileups() {
		List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(cellBCFile);
		IntervalList snpIntervals = IntervalList.fromFile(snpIntervalsFile);
		
		// all SNPs equal quality.  No SNPs overlap on UMIs so all pileups are seen.
		Map<Interval, Double> meanGenotypeQuality = new HashMap<Interval, Double>();		
		for (Interval i: snpIntervals) {
			meanGenotypeQuality.put(i, 30d);
		}
		
		SNPUMIBasePileupIterator sbpi = new SNPUMIBasePileupIterator(
				new SamHeaderAndIterator(smallBAMFile), snpIntervals, GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG,
				LOCUS_FUNCTION_LIST, STRAND_STRATEGY, cellBarcodeTag,
				molBCTag, snpTag, functionTag, readMQ, assignReadsToAllGenes, cellBarcodes, meanGenotypeQuality, SortOrder.SNP_GENE);

		int expectedPileups=8;
		int count=0;
		while (sbpi.hasNext()) {
			SNPUMIBasePileup p = sbpi.next();
			checkAnswer(p);
			count++;
		}
		Assert.assertEquals(expectedPileups, count);
	}
	
	
	// for reads that have more than 1 SNP in the data, filter the SNPs tagging the read to a single SNP.
	// This means the read only contributes to a single pileup.
	@Test(enabled=true)
	public void testDuplicateReadsWithMultipleSNPs() {
		List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(cellBCFile);
		IntervalList snpIntervals = IntervalList.fromFile(snpIntervalsFile);
		// add an interval one base later.  This intersects a subset of reads.  
		snpIntervals.add(new Interval("HUMAN_1", 76227021, 76227021, false, "reject"));
		
		
		// if we mark the 3nd snp as lower quality we have 8 pileups, with the snp named HUMAN_1:76227022
		Map<Interval, Double> genotypeQuality = new HashMap<>();
		genotypeQuality.put(snpIntervals.getIntervals().get(0), 100d);
		genotypeQuality.put(snpIntervals.getIntervals().get(1), 50d);
		genotypeQuality.put(snpIntervals.getIntervals().get(2), 10d);
		testSNPQualities (snpIntervals, cellBarcodes, genotypeQuality, 8, "HUMAN_1:76227022");
		
		// if we mark the 1st snp as lower quality we have 8 pileups. The snp named "reject" is in all the pileups	
		genotypeQuality = new HashMap<>();
		genotypeQuality.put(snpIntervals.getIntervals().get(0), 10d);		
		genotypeQuality.put(snpIntervals.getIntervals().get(2), 100d);
		testSNPQualities (snpIntervals, cellBarcodes, genotypeQuality, 8, "reject");
				
	}

	private void testSNPQualities (IntervalList snpIntervals, List<String> cellBarcodes, Map<Interval, Double> genotypeQuality, int expectedPileUps, String expectedSnpNamePrefix) {
		SNPUMIBasePileupIterator sbpi = new SNPUMIBasePileupIterator(
				new SamHeaderAndIterator(smallBAMFile), snpIntervals, GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG,
				LOCUS_FUNCTION_LIST, STRAND_STRATEGY, cellBarcodeTag,
				molBCTag, snpTag, functionTag, readMQ, assignReadsToAllGenes, cellBarcodes, genotypeQuality, SortOrder.SNP_GENE);
		
		int count=0;
		while (sbpi.hasNext()) {
			SNPUMIBasePileup p = sbpi.next();
			if (expectedSnpNamePrefix!=null) 
				Assert.assertTrue(p.getSnpID().contains(expectedSnpNamePrefix));
			System.out.println (p);
			count++;
		}
		System.out.println ();
		Assert.assertEquals(expectedPileUps, count);		
	}

	private void checkAnswer (final SNPUMIBasePileup p) {
		// G37, G37, G27
		if (p.getCell().equals("ATCAGGGACAGA") && p.getMolecularBarcode().equals("CGGGGCTC")) {
			char [] eb = {'G', 'G', 'G'};
			byte [] eq = {37,27,37};
			testPileup(p, eb, eq);
		}

		if (p.getCell().equals("ATCAGGGACAGA") && p.getMolecularBarcode().equals("CTGGGCTC")) {
			char [] eb = {'G'};
			byte [] eq = {37};
			testPileup(p, eb, eq);
		}

		if (p.getCell().equals("ATCAGGGACAGA") && p.getMolecularBarcode().equals("GGAATGTG")) {
			char [] eb = {'A'};
			byte [] eq = {8};
			testPileup(p, eb, eq);
		}

		if (p.getCell().equals("ATCAGGGACAGA") && p.getMolecularBarcode().equals("GTGCCGTG")) {
			char [] eb = {'A'};
			byte [] eq = {8};
			testPileup(p, eb, eq);
		}

		if (p.getCell().equals("ATCAGGGACAGA") && p.getMolecularBarcode().equals("TCGCAGAC")) {
			char [] eb = {'G'};
			byte [] eq = {37};
			testPileup(p, eb, eq);
		}

		if (p.getCell().equals("TTGCCTTACGCG") && p.getMolecularBarcode().equals("CAGCGGTG")) {
			char [] eb = {'G', 'G'};
			byte [] eq = {37,32};
			testPileup(p, eb, eq);
		}

		if (p.getCell().equals("TTGCCTTACGCG") && p.getMolecularBarcode().equals("GGTGTCAT")) {
			char [] eb = {'A', 'A', 'A'};
			byte [] eq = {27,32,13};
			testPileup(p, eb, eq);
		}

		if (p.getCell().equals("TTGCCTTACGCG") && p.getMolecularBarcode().equals("TGGTGGGG")) {
			char [] eb = {'A'};
			byte [] eq = {37};
			testPileup(p, eb, eq);
		}


	}

	private void testPileup (final SNPUMIBasePileup p, final char [] expectedBases, final byte [] expectedQuals) {
		List<Character> bases = p.getBasesAsCharacters();

		List<Byte> quals = p.getQualities();

		char[] basesArray = getAsCharArray(bases);
		byte[] qualsArray = getAsByteArray(quals);
		Assert.assertEquals(expectedBases, basesArray);
		Assert.assertEquals(expectedQuals, qualsArray);
	}

	private char [] getAsCharArray (final List<Character> l) {
		char [] result = new char[l.size()];
		for (int i=0; i<result.length; i++)
			result[i]=l.get(i);
		return (result);
	}

	private byte [] getAsByteArray (final List<Byte> l) {
		byte [] result = new byte[l.size()];
		for (int i=0; i<result.length; i++)
			result[i]=l.get(i);
		return (result);
	}



}
