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

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPUMIBasePileupIterator;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SortOrder;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IntervalList;
import picard.annotation.LocusFunction;

public class SampleGenotypeProbabilitiesIteratorTest {

	// See smallTest_snpUMIPileUp.summary.txt for a listing of the
	// bases/qualities/cell/molecular barcodes for all 15 reads
	
	private final static String rootDir="testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts";
	
	private final File INPUT_BAM = new File(rootDir + "/smallTest_snpUMIPileUp_retagged.sam");
	private final File snpIntervalsFile = new File(rootDir + "/hek_cells_2snps.intervals");


	private final String cellBarcodeTag = "XC";
	private final String molBCTag = "XM";
	private final String snpTag = "YS";
	private final int readMQ = 10;
	private final boolean assignReadsToAllGenes = true;

	private String GENE_NAME_TAG="gn";
	private String GENE_STRAND_TAG="gs";
	private String GENE_FUNCTION_TAG="gf";
	private StrandStrategy STRAND_STRATEGY = StrandStrategy.SENSE;
	private List<LocusFunction> LOCUS_FUNCTION_LIST=new ArrayList<LocusFunction>(Arrays.asList(LocusFunction.CODING, LocusFunction.UTR));

	@Test
	public void testIterED1 () {
		int editDistance=1;
		List<String> barcodes = new ArrayList<>();
		barcodes.add("ATCAGGGACAGA");
		SampleGenotypeProbabilitiesIterator iter = prepIterFromFile(barcodes, editDistance);
		// only 1 entry.
		int count=0;
		while (iter.hasNext()) {
			SampleGenotypeProbabilities p = iter.next();
			int numPileups = p.getBackingPileups().size();
			Assert.assertEquals(numPileups, 4);
			count++;
		}
		Assert.assertEquals(1, count);
	}

	@Test
	public void testIterED0 () {
		int editDistance=0;
		List<String> barcodes = new ArrayList<>();
		barcodes.add("ATCAGGGACAGA");
		SampleGenotypeProbabilitiesIterator iter = prepIterFromFile(barcodes, editDistance);
		// only 1 entry.
		int count=0;
		while (iter.hasNext()) {
			SampleGenotypeProbabilities p = iter.next();
			int numPileups = p.getBackingPileups().size();
			Assert.assertEquals(numPileups, 5);
			count++;
		}
		Assert.assertEquals(1, count);
	}

	@Test
	public void testTwoCells () {
		int editDistance=0;
		List<String> barcodes = new ArrayList<>();
		barcodes.add("ATCAGGGACAGA");
		barcodes.add("TTGCCTTACGCG");
		SampleGenotypeProbabilitiesIterator iter = prepIterFromFile(barcodes, editDistance);
		// only 1 entry.
		int count=0;
		while (iter.hasNext()) {
			SampleGenotypeProbabilities p = iter.next();
			testTwoCellResults(p);
			count++;
		}
		Assert.assertEquals(2, count);
	}

	private void testTwoCellResults (final SampleGenotypeProbabilities p) {
		if (p.getCell().equals("ATCAGGGACAGA")) {
			int numPileups = p.getBackingPileups().size();
			Assert.assertEquals(numPileups, 5);
		}
		if (p.getCell().equals("TTGCCTTACGCG")) {
			int numPileups = p.getBackingPileups().size();
			Assert.assertEquals(numPileups, 3);
		}
	}


	private SampleGenotypeProbabilitiesIterator prepIterFromFile (final List<String> barcodes, final int editDistance) {
		IntervalList snpIntervals = IntervalList.fromFile(snpIntervalsFile);

		SNPUMIBasePileupIterator sbpi = new SNPUMIBasePileupIterator(
				new SamHeaderAndIterator(this.INPUT_BAM), snpIntervals, GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG,
				LOCUS_FUNCTION_LIST, STRAND_STRATEGY, cellBarcodeTag,
				molBCTag, snpTag, "XF", readMQ, assignReadsToAllGenes, barcodes, null, SortOrder.SNP_CELL);

		boolean flag = sbpi.hasNext();
		final SAMSequenceDictionary dict = snpIntervals.getHeader().getSequenceDictionary();
		SampleGenotypeProbabilitiesIterator result = new SampleGenotypeProbabilitiesIterator(sbpi, dict, editDistance);
		return result;
	}
}
