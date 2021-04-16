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
import java.util.Collection;
import java.util.List;

import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import org.junit.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.util.IntervalList;
import picard.annotation.LocusFunction;

public class MultiCellDigitalAlleleCountsIteratorTest {

	// See smallTest_snpUMIPileUp.summary.txt for a listing of the bases/qualities/cell/molecular barcodes for all 15 reads
		private final File smallBAMFile = new File(
				"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/smallTest_snpUMIPileUp_retagged.sam");
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

		private File largeBAMFile = new File(
				"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/hek_5_cell_2_snp_testdata_retagged.bam");

		@Test(enabled=true)
		public void testGatherMultiCellDACs() {
			List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(cellBCFile);
			IntervalList snpIntervals = IntervalList.fromFile(snpIntervalsFile);

			int baseQualityThreshold=10;
			SNPUMIBasePileupIterator sbpi = new SNPUMIBasePileupIterator(
					new SamHeaderAndIterator(smallBAMFile), snpIntervals, GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG,
					LOCUS_FUNCTION_LIST, STRAND_STRATEGY, cellBarcodeTag,
					molBCTag, snpTag, functionTag, readMQ, assignReadsToAllGenes, cellBarcodes, null, SortOrder.SNP_GENE);

			MultiCellDigitalAlleleCountsIterator multiIter = new MultiCellDigitalAlleleCountsIterator(new DigitalAlleleCountsIterator(sbpi, baseQualityThreshold));

			int counter=0;
			while (multiIter.hasNext()) {
				MultiCellDigitalAlleleCounts mcdac = multiIter.next();
				for (String cellID: mcdac.getCells()) {
					DigitalAlleleCounts dac = mcdac.getDigitalAlleleCounts(cellID);
					checkAnswerBQ10(dac);
				}
				counter++;
			}
			// assert that you saw 1 gene/snp MultiCellDAC that has 0 or more cells.
			Assert.assertEquals(1, counter);
			multiIter.close();
		}

		/**
		 * For the smallTest_snpUMIPileUp data set, there are 2 DACs.
		 * Base quality will change this answer.
		 * Only test the reads and umis observed.  The derived values from these counts should be tested by DigitalAlleleCountsTest
		 * @param dac
		 */
		private void checkAnswerBQ10 (final DigitalAlleleCounts dac) {
			String cell = dac.getCell();
			if (cell.equals("ATCAGGGACAGA")) {
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

		/**
		 * Using BQ=0 here to make data similar to the ad-hoc analysis that started this.
		 * This data set assumes no UMI collapse.
		 * 5 cells, 2 SNPs + meta analysis.
		 */
		@Test(enabled=true)
		public void bigDataTest () {
			List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(cellBCFile);
			IntervalList snpIntervals = IntervalList.fromFile(snpIntervalsFile);
			int baseQualityThreshold=0;

			SNPUMIBasePileupIterator sbpi = new SNPUMIBasePileupIterator(
					new SamHeaderAndIterator(largeBAMFile), snpIntervals, GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG,
					LOCUS_FUNCTION_LIST, STRAND_STRATEGY, cellBarcodeTag,
					molBCTag, snpTag, functionTag, readMQ, assignReadsToAllGenes, cellBarcodes, null, SortOrder.SNP_GENE);

			MultiCellDigitalAlleleCountsIterator multiIter = new MultiCellDigitalAlleleCountsIterator(new DigitalAlleleCountsIterator(sbpi, baseQualityThreshold));

			int counter=0;
			while (multiIter.hasNext()) {
				MultiCellDigitalAlleleCounts mcdac = multiIter.next();
				for (String cellID: mcdac.getCells()) {
					DigitalAlleleCounts dac = mcdac.getDigitalAlleleCounts(cellID);
					Collection<String> umis = dac.umis();
					int numUMIs=dac.umis().size();
					checkAnswerLargeData(dac);
				}
				DigitalAlleleCounts meta = mcdac.getMetaAnalysis();
				checkAnswerLargeData(meta);
				counter++;
			}
			// assert that you saw 1 gene/snp MultiCellDAC that has 0 or more cells.
			Assert.assertEquals(2, counter);
			multiIter.close();
		}

		/**
		 * Using BQ=0 here to make data similar to the ad-hoc analysis that started this.
		 * This data set has an ED=0 collapse, which should be the equivalent of the bigDataTest results.
		 * 5 cells, 2 SNPs + meta analysis.
		 */
		@Test(enabled=true)
		public void bigDataTestED0 () {
			List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(cellBCFile);
			IntervalList snpIntervals = IntervalList.fromFile(snpIntervalsFile);

			int baseQualityThreshold=0;

			SNPUMIBasePileupIterator sbpi = new SNPUMIBasePileupIterator(
					new SamHeaderAndIterator(largeBAMFile), snpIntervals, GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG,
					LOCUS_FUNCTION_LIST, STRAND_STRATEGY, cellBarcodeTag,
					molBCTag, snpTag, functionTag, readMQ, assignReadsToAllGenes, cellBarcodes, null, SortOrder.SNP_GENE);

			MultiCellDigitalAlleleCountsIterator multiIter = new MultiCellDigitalAlleleCountsIterator(new DigitalAlleleCountsIterator(sbpi, baseQualityThreshold));

			int counter=0;
			while (multiIter.hasNext()) {
				MultiCellDigitalAlleleCounts mcdac = multiIter.next();
				mcdac.collapseDACs(0);
				for (String cellID: mcdac.getCells()) {
					DigitalAlleleCounts dac = mcdac.getDigitalAlleleCounts(cellID);
					Collection<String> umis = dac.umis();
					int numUMIs=dac.umis().size();
					checkAnswerLargeData(dac);
				}
				DigitalAlleleCounts meta = mcdac.getMetaAnalysis();
				checkAnswerLargeData(meta);
				counter++;
			}
			// assert that you saw 1 gene/snp MultiCellDAC that has 0 or more cells.
			Assert.assertEquals(2, counter);
			multiIter.close();
		}


		private void checkAnswerLargeData(final DigitalAlleleCounts dac) {
			if (dac.getSnpInterval().getName().equals("HUMAN_1:150199123") && dac.getCell().equals("AGGGAAAATTGA"))
				checkAnswerSingle(dac, 'G', 'T', 43, 14, 15, 4);

			if (dac.getSnpInterval().getName().equals("HUMAN_1:150199123") && dac.getCell().equals("ATCAGGGACAGA"))
				checkAnswerSingle(dac, 'G', 'T', 62, 46, 17, 22);

			if (dac.getSnpInterval().getName().equals("HUMAN_1:150199123") && dac.getCell().equals("TACAATTAAGGC"))
				checkAnswerSingle(dac, 'G', 'T', 43, 2, 8, 2);

			if (dac.getSnpInterval().getName().equals("HUMAN_1:150199123") && dac.getCell().equals("TGGCGAAGAGAT"))
				checkAnswerSingle(dac, 'G', 'T', 12, 7, 5, 3);

			if (dac.getSnpInterval().getName().equals("HUMAN_1:150199123") && dac.getCell().equals("TTGCCTTACGCG"))
				checkAnswerSingle(dac, 'G', 'T', 52, 5, 11, 1);

			if (dac.getSnpInterval().getName().equals("HUMAN_1:150199123") && dac.getCell().equals("ALL_CELLS"))
				checkAnswerSingle(dac, 'G', 'T', 212,74,56,32);


			if (dac.getSnpInterval().getName().equals("HUMAN_1:76227022") && dac.getCell().equals("AGGGAAAATTGA"))
				checkAnswerSingle(dac, 'A', 'G', 53, 40, 8, 8);

			if (dac.getSnpInterval().getName().equals("HUMAN_1:76227022") && dac.getCell().equals("ATCAGGGACAGA"))
				checkAnswerSingle(dac, 'A', 'G', 85, 72, 17, 14);

			if (dac.getSnpInterval().getName().equals("HUMAN_1:76227022") && dac.getCell().equals("TACAATTAAGGC"))
				checkAnswerSingle(dac, 'A', 'G', 28, 5, 5, 3);

			if (dac.getSnpInterval().getName().equals("HUMAN_1:76227022") && dac.getCell().equals("TGGCGAAGAGAT"))
				checkAnswerSingle(dac, 'A', 'G', 24, 2, 5, 2);

			if (dac.getSnpInterval().getName().equals("HUMAN_1:76227022") && dac.getCell().equals("TTGCCTTACGCG"))
				checkAnswerSingle(dac, 'A', 'G', 61,65,13,9);

			if (dac.getSnpInterval().getName().equals("HUMAN_1:76227022") && dac.getCell().equals("ALL_CELLS"))
				checkAnswerSingle(dac, 'A', 'G', 251,184,48,36);
		}



		/**
		 * Using BQ=0 here to make data similar to the ad-hoc analysis that started this.
		 * This data set has an ED=0 collapse, which should be the equivalent of the bigDataTest results.
		 * 5 cells, 2 SNPs + meta analysis.
		 */
		@Test(enabled=true)
		public void bigDataTestED1 () {
			List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(cellBCFile);
			IntervalList snpIntervals = IntervalList.fromFile(snpIntervalsFile);
			int baseQualityThreshold=0;
			int editDistance=1;
			SNPUMIBasePileupIterator sbpi = new SNPUMIBasePileupIterator(
					new SamHeaderAndIterator(largeBAMFile), snpIntervals, GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG,
					LOCUS_FUNCTION_LIST, STRAND_STRATEGY, cellBarcodeTag,
					molBCTag, snpTag, functionTag, readMQ, assignReadsToAllGenes, cellBarcodes, null, SortOrder.SNP_GENE);

			MultiCellDigitalAlleleCountsIterator multiIter = new MultiCellDigitalAlleleCountsIterator(new DigitalAlleleCountsIterator(sbpi, baseQualityThreshold));

			int counter=0;
			while (multiIter.hasNext()) {
				MultiCellDigitalAlleleCounts mcdac = multiIter.next();
				mcdac.collapseDACs(editDistance);
				for (String cellID: mcdac.getCells()) {
					DigitalAlleleCounts dac = mcdac.getDigitalAlleleCounts(cellID);
					Collection<String> umis = dac.umis();
					int numUMIs=dac.umis().size();
					checkAnswerLargeDataED1(dac);
				}
				DigitalAlleleCounts meta = mcdac.getMetaAnalysis();
				checkAnswerLargeDataED1(meta);
				counter++;
			}
			// assert that you saw 1 gene/snp MultiCellDAC that has 0 or more cells.
			Assert.assertEquals(2, counter);
			multiIter.close();
		}

		public void checkAnswerLargeDataED1(final DigitalAlleleCounts dac) {
			//[{T=14, G=43}] UMI counts [{T=4, G=14}]
			if (dac.getSnpInterval().getName().equals("HUMAN_1:150199123") && dac.getCell().equals("AGGGAAAATTGA"))
				checkAnswerSingle(dac, 'G', 'T', 43, 14, 14, 4);
			//[{T=46, G=62}] UMI counts [{T=19, G=16}]
			if (dac.getSnpInterval().getName().equals("HUMAN_1:150199123") && dac.getCell().equals("ATCAGGGACAGA"))
				checkAnswerSingle(dac, 'G', 'T', 62, 46, 16, 19);
			// [{T=2, G=43}] UMI counts [{T=2, G=8}]
			if (dac.getSnpInterval().getName().equals("HUMAN_1:150199123") && dac.getCell().equals("TACAATTAAGGC"))
				checkAnswerSingle(dac, 'G', 'T', 43, 2, 8, 2);
			//  [{T=7, G=12}] UMI counts [{T=3, G=5}]
			if (dac.getSnpInterval().getName().equals("HUMAN_1:150199123") && dac.getCell().equals("TGGCGAAGAGAT"))
				checkAnswerSingle(dac, 'G', 'T', 12, 7, 5, 3);
			// [{T=5, G=52}] UMI counts [{T=1, G=11}]
			if (dac.getSnpInterval().getName().equals("HUMAN_1:150199123") && dac.getCell().equals("TTGCCTTACGCG"))
				checkAnswerSingle(dac, 'G', 'T', 52, 5, 11, 1);

			// [{T=74, G=212}] UMI counts [{T=29, G=54}]
			if (dac.getSnpInterval().getName().equals("HUMAN_1:150199123") && dac.getCell().equals("ALL_CELLS"))
				checkAnswerSingle(dac, 'G', 'T', 212,74,54,29);

			// [{G=40, A=53}] UMI counts [{G=6, A=8}]
			if (dac.getSnpInterval().getName().equals("HUMAN_1:76227022") && dac.getCell().equals("AGGGAAAATTGA"))
				checkAnswerSingle(dac, 'A', 'G', 53, 40, 8, 6);
			// [{T=1, G=72, A=85}] UMI counts [{G=13, A=14}]
			if (dac.getSnpInterval().getName().equals("HUMAN_1:76227022") && dac.getCell().equals("ATCAGGGACAGA"))
				checkAnswerSingle(dac, 'A', 'G', 85, 72, 14, 13);
			// [{G=5, A=28}] UMI counts [{G=3, A=4}]
			if (dac.getSnpInterval().getName().equals("HUMAN_1:76227022") && dac.getCell().equals("TACAATTAAGGC"))
				checkAnswerSingle(dac, 'A', 'G', 28, 5, 4, 3);
			// [{G=2, A=24}] UMI counts [{G=2, A=5}]
			if (dac.getSnpInterval().getName().equals("HUMAN_1:76227022") && dac.getCell().equals("TGGCGAAGAGAT"))
				checkAnswerSingle(dac, 'A', 'G', 24, 2, 5, 2);
			// [{G=65, A=61}] UMI counts [{G=9, A=10}]
			if (dac.getSnpInterval().getName().equals("HUMAN_1:76227022") && dac.getCell().equals("TTGCCTTACGCG"))
				checkAnswerSingle(dac, 'A', 'G', 61,65,10,9);
			// [{T=1, G=184, A=251}] UMI counts [{G=33, A=41}]
			if (dac.getSnpInterval().getName().equals("HUMAN_1:76227022") && dac.getCell().equals("ALL_CELLS"))
				checkAnswerSingle(dac, 'A', 'G', 251,184,41,33);


		}

		private void checkAnswerSingle (final DigitalAlleleCounts dac, final char alleleOne, final char alleleTwo, final int readCountExpectedA1, final int readCountExpectedA2, final int umiCountExpectedA1, final int umiCountExpectedA2) {
			ObjectCounter<Character> reads = dac.getReadCounts();
			int readCountObservedA1=reads.getCountForKey(alleleOne);
			int readCountObservedA2=reads.getCountForKey(alleleTwo);
			Assert.assertEquals(readCountExpectedA1, readCountObservedA1);
			Assert.assertEquals(readCountExpectedA2, readCountObservedA2);
			ObjectCounter<Character> umis = dac.getUMIAlleleCount();
			int umiCountObserveddA1=umis.getCountForKey(alleleOne);
			int umiCountObserveddA2=umis.getCountForKey(alleleTwo);
			Assert.assertEquals(umiCountExpectedA1, umiCountObserveddA1);
			Assert.assertEquals(umiCountExpectedA2, umiCountObserveddA2);
		}


}
