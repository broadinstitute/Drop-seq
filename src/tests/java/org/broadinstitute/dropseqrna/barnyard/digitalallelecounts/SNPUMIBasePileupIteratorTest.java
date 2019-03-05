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

import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPUMIBasePileup;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPUMIBasePileupIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import org.junit.Assert;
import org.testng.annotations.Test;
import picard.annotation.LocusFunction;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SNPUMIBasePileupIteratorTest {

	// See smallTest_snpUMIPileUp.summary.txt for a listing of the bases/qualities/cell/molecular barcodes for all 15 reads
	private final File smallBAMFile = new File(
			"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/smallTest_snpUMIPileUp.sam");
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

		SNPUMIBasePileupIterator sbpi = new SNPUMIBasePileupIterator(
				smallBAMFile, snpIntervals, GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG,
				LOCUS_FUNCTION_LIST, STRAND_STRATEGY, cellBarcodeTag,
				molBCTag, snpTag, functionTag, readMQ, assignReadsToAllGenes, cellBarcodes);

		while (sbpi.hasNext()) {
			SNPUMIBasePileup p = sbpi.next();
			checkAnswer(p);
		}

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
		Assert.assertArrayEquals(expectedBases, basesArray);
		Assert.assertArrayEquals(expectedQuals, qualsArray);
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
