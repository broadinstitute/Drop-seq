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

import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalDataProcessorStrategy;
import org.broadinstitute.dropseqrna.barnyard.GeneFunctionCommandLineBase;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import picard.annotation.LocusFunction;

public class DigitalAlleleCountsIteratorTest {

	// See smallTest_snpUMIPileUp.summary.txt for a listing of the bases/qualities/cell/molecular barcodes for all 15 reads
	private final File smallBAMFile = new File(
			"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/smallTest_snpUMIPileUp_retagged.sam");
	
	private final File smallBAMFile2 = new File (
			"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/smallTest_snpUMIPileUp_retagged_bestGeneTest.sam");
	
	private final File cellBCFile = new File(
			"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/hek_cells_cell_barcodes.txt");
	private final File snpIntervalsFile = new File(
			"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/hek_cells_2snps.intervals");
	private final String cellBarcodeTag = "XC";
	private final String molBCTag="XM";
	private String GENE_NAME_TAG="gn";
	private String GENE_STRAND_TAG="gs";
	private String GENE_FUNCTION_TAG="gf";
	private StrandStrategy STRAND_STRATEGY = StrandStrategy.SENSE;
	private List<LocusFunction> LOCUS_FUNCTION_LIST=new ArrayList<LocusFunction>(Arrays.asList(LocusFunction.CODING, LocusFunction.UTR));private final String snpTag = "YS";
    private final String functionTag="XF";
	private final int readMQ = 10;
	private final boolean assignReadsToAllGenes = true;

	
	// An adapter now that meanGenotypeQuality is required.
	private Map<Interval, Double> getAverageGQ (IntervalList snpIntervals) {
		// set all genotype qualities to missing.
		Map<Interval,Double>meanGenotypeQuality =new HashMap<>();
		for (Interval i: snpIntervals.getIntervals())
			meanGenotypeQuality.put(i,-1d);
		return meanGenotypeQuality;
	}
	
	@Test
	public void testGatherDACs() {
		List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(cellBCFile);
		IntervalList snpIntervals = IntervalList.fromFile(snpIntervalsFile);
		int baseQualityThreshold=10;
		Map<Interval, Double> meanGQ= getAverageGQ(snpIntervals);

		SNPUMIBasePileupIterator sbpi = new SNPUMIBasePileupIterator(
				new SamHeaderAndIterator(smallBAMFile), snpIntervals, GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG,
				LOCUS_FUNCTION_LIST,STRAND_STRATEGY, FunctionalDataProcessorStrategy.DROPSEQ, cellBarcodeTag, molBCTag, snpTag,
				GeneFunctionCommandLineBase.DEFAULT_FUNCTION_TAG, readMQ, assignReadsToAllGenes,
				cellBarcodes, meanGQ, SortOrder.SNP_GENE);

		DigitalAlleleCountsIterator daci = new DigitalAlleleCountsIterator(sbpi, baseQualityThreshold);
		int counter=0;
		while (daci.hasNext()) {
			DigitalAlleleCounts dac = daci.next();
			checkAnswerBQ10(dac);
			counter++;
		}
		// assert that you saw 2 dacs, for the 2 cells.
		Assert.assertEquals(2, counter);
	}
	
	
	
	@Test
	// this has additional reads on a new gene called "FAKE"
	public void testGatherDACs2() {
		List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(cellBCFile);
		IntervalList snpIntervals = IntervalList.fromFile(snpIntervalsFile);
		int baseQualityThreshold=10;
		
		Map<Interval, Double> meanGQ= getAverageGQ(snpIntervals);
		
		SNPUMIBasePileupIterator sbpi = new SNPUMIBasePileupIterator(
				new SamHeaderAndIterator(smallBAMFile2), snpIntervals, GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG,
				LOCUS_FUNCTION_LIST,STRAND_STRATEGY, FunctionalDataProcessorStrategy.DROPSEQ, cellBarcodeTag, molBCTag, snpTag,
				GeneFunctionCommandLineBase.DEFAULT_FUNCTION_TAG, readMQ, assignReadsToAllGenes,
				cellBarcodes, meanGQ, SortOrder.SNP_GENE);

		DigitalAlleleCountsIterator daci = new DigitalAlleleCountsIterator(sbpi, baseQualityThreshold);
		int counter=0;
		while (daci.hasNext()) {
			DigitalAlleleCounts dac = daci.next();
			checkAnswerBQ10(dac);
			counter++;
		}
		// assert that you saw 3 dacs, for the 2 cells.
		Assert.assertEquals(3, counter);
	}
	
	
	/**
	 * For the smallTest_snpUMIPileUp data set, there are 2 DACs.
	 * Base quality will change this answer.
	 * Only test the reads and umis observed.  The derived values from these counts should be tested by DigitalAlleleCountsTest
	 * @param dac
	 */
	private void checkAnswerBQ10 (final DigitalAlleleCounts dac) {
		String cell = dac.getCell();
		if (cell.equals("ATCAGGGACAGA") && dac.getGene().equals("ACADM")) {
			String snpID = dac.getSnpInterval().getName();
			Assert.assertEquals("HUMAN_1:76227022", snpID);
			String gene = dac.getGene();
			Assert.assertEquals("ACADM", gene);

			// test read counts
			ObjectCounter<Character> readCounts= dac.getReadCounts();
			int readG = readCounts.getCountForKey('G');
			int readA = readCounts.getCountForKey('A');
			Assert.assertEquals(5, readG);
			Assert.assertEquals(0, readA);

			// test umi counts
			ObjectCounter<Character> umiCounts= dac.getUMIAlleleCount();
			int umiG = umiCounts.getCountForKey('G');
			int umiA = umiCounts.getCountForKey('A');
			Assert.assertEquals(3, umiG);
			Assert.assertEquals(0, umiA);
		}
		
		if (cell.equals("ATCAGGGACAGA") && dac.getGene().equals("FAKE")) {
			String snpID = dac.getSnpInterval().getName();
			Assert.assertEquals("HUMAN_1:76227022", snpID);
			String gene = dac.getGene();
			Assert.assertEquals("FAKE", gene);

			// test read counts
			ObjectCounter<Character> readCounts= dac.getReadCounts();
			int readG = readCounts.getCountForKey('G');
			int readA = readCounts.getCountForKey('A');
			Assert.assertEquals(4, readG);
			Assert.assertEquals(0, readA);

			// test umi counts
			ObjectCounter<Character> umiCounts= dac.getUMIAlleleCount();
			int umiG = umiCounts.getCountForKey('G');
			int umiA = umiCounts.getCountForKey('A');
			Assert.assertEquals(2, umiG);
			Assert.assertEquals(0, umiA);
		}

		if (cell.equals("TTGCCTTACGCG")) {
			String snpID = dac.getSnpInterval().getName();
			Assert.assertEquals("HUMAN_1:76227022", snpID);
			String gene = dac.getGene();
			Assert.assertEquals("ACADM", gene);

			// test read counts
			ObjectCounter<Character> readCounts= dac.getReadCounts();
			int readG = readCounts.getCountForKey('G');
			int readA = readCounts.getCountForKey('A');
			Assert.assertEquals(2, readG);
			Assert.assertEquals(4, readA);

			// test umi counts
			ObjectCounter<Character> umiCounts= dac.getUMIAlleleCount();
			int umiG = umiCounts.getCountForKey('G');
			int umiA = umiCounts.getCountForKey('A');
			Assert.assertEquals(1, umiG);
			Assert.assertEquals(2, umiA);
		}


	}


}
