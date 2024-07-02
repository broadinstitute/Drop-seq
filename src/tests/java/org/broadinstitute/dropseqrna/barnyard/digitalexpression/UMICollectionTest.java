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

import htsjdk.samtools.util.IterableAdapter;
import org.broadinstitute.dropseqrna.barnyard.DigitalExpressionTest;
import org.broadinstitute.dropseqrna.barnyard.DigitalExpressionTestUtil;
import org.broadinstitute.dropseqrna.barnyard.GeneFunctionCommandLineBase;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.sam.util.Pair;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.File;
import java.util.*;


public class UMICollectionTest {

	File umiCollectionFole = new File("testdata/org/broadinstitute/transcriptome/barnyard/digitalexpression/UMICollectionFile.txt.gz");


	@Test()
	public void test1 () {

		UMICollection c = new UMICollection("C", "G");
		c.incrementMolecularBarcodeCount("FOO");
		c.incrementMolecularBarcodeCount("BAR");
		c.incrementMolecularBarcodeCount("BAH");

		int dge = c.getDigitalExpression(1, 1, false);
		Assert.assertEquals(dge, 2);

		int count = c.getDigitalExpression(1, 1, true);
		Assert.assertEquals(count, 3);

		int dge2Min=c.getDigitalExpression(2, 1, false);
		Assert.assertEquals(dge2Min, 1);

		String r = c.toString();
		String expected ="[C] + [G] {BAR=1, FOO=1, BAH=1}";
		Assert.assertNotNull(r);
		Assert.assertEquals(r, expected);

		ObjectCounter <String> counts = c.getMolecularBarcodeCounts();
		Assert.assertNotNull(counts);
		Assert.assertEquals(counts.toString(), "{BAR=1, FOO=1, BAH=1}");

	}

	@Test()
	public void testFileParse () {
		Collection<UMICollection> result=  UMICollection.parseUMICollectionFile(this.umiCollectionFole);

		Assert.assertNotNull(result);
		Assert.assertEquals(result.size(), 12);

	}

	@Test()
	public void testFilter () {
		// model a set that will get filtered.
		UMICollection c = new UMICollection("C", "G");
		c.incrementMolecularBarcodeCount("FOO", 1000);
		c.incrementMolecularBarcodeCount("BAR", 1000);
		c.incrementMolecularBarcodeCount("BAH", 1);
		c.filterByUMIFrequency(0.01);
		Assert.assertEquals(c.getDigitalExpression(0, 0, false), 2);

		// model a set that will NOT get filtered.
		UMICollection c2 = new UMICollection("C", "G");
		c2.incrementMolecularBarcodeCount("FOO", 10);
		c2.incrementMolecularBarcodeCount("BAR", 100);
		c2.incrementMolecularBarcodeCount("BAH", 1);
		c2.filterByUMIFrequency(0.01);
		Assert.assertEquals(c2.getDigitalExpression(0, 0, false), 3);


	}

    private void testDownsampleRate(double downsampleRate, double numReads, double confidenceRate) {
        /* Suppose each read of a UMI has a 0.5 chance of getting dropped
         * This means that we expect 50 in 100 counts to be dropped
         * with some standard deviation (say 5%?)
         */
        UMICollection c = new UMICollection("C", "G");
        Random random = new Random(1L);
        int nr=(int)(numReads/4);
        int nr2 = (int)(numReads - ((double)nr*3));
        c.incrementMolecularBarcodeCount("FOO", nr);
        c.incrementMolecularBarcodeCount("BAR", nr);
        c.incrementMolecularBarcodeCount("FOOBAR", nr);
        c.incrementMolecularBarcodeCount("BARFOO", nr2);
        ObjectCounter<String> counts = c.getDownsampledMolecularBarcodeCounts(downsampleRate, random);
        Assert.assertEquals(counts.getTotalCount(), numReads*downsampleRate, numReads*confidenceRate);
    }

    @Test()
    public void testDownsample() {
        /* Make new UMICollection object
         * Initialize counts
         * Test downsampling at different rates
         */
        //testDownsampleRate(downsampleRate, 100.0, 0.05);
        double numReads=100.0;
        double confidenceRate=0.05;
        testDownsampleRate( 0.0, numReads, confidenceRate);
        testDownsampleRate( 1.0, numReads, confidenceRate);
        testDownsampleRate( 0.5, numReads, confidenceRate);
        testDownsampleRate(0.75, numReads, confidenceRate);
        testDownsampleRate(0.001, 1000000.0, 0.0005);
    }

	// Do UMI collapse and confirm that the UMI counts are the same as DigitalExpressionTest.doWork()
	@Test
	public void testEditDistanceCollapse() {
		final UMIIterator umiIterator = new UMIIterator.UMIIteratorBuilder(SamFileMergeUtil.mergeInputs(Collections.singletonList(DigitalExpressionTestUtil.IN_FILE), false),
				DigitalExpressionTest.GENE_NAME_TAG,
				DigitalExpressionTest.GENE_STRAND_TAG, DigitalExpressionTest.GENE_FUNCTION_TAG, DigitalExpressionTest.STRAND_STRATEGY,
				DigitalExpressionTest.LOCUS_FUNCTION_LIST, GeneFunctionCommandLineBase.DEFAULT_FUNCTIONAL_STRATEGY,
				"XC", DigitalExpressionTest.MOLECULAR_BARCODE_TAG, DigitalExpressionTest.READ_MQ).
				setCellBarcodes(Arrays.asList(DigitalExpressionTestUtil.barcodes)).retainReads(true).build();

		final Map<Pair<String, String>, Integer> umiCounts = new HashMap<>();
		for (final UMICollection c : new IterableAdapter<UMICollection>(umiIterator)) {
			c.collapseThisByEditDistance(1, DigitalExpressionTest.MOLECULAR_BARCODE_TAG);
			final Pair<String, String> key = new Pair<>(c.getCellBarcode(), c.getGeneName());
			umiCounts.put(key, c.getMolecularBarcodeCounts().getSize());
		}
		final Map<Pair<String, String>, Integer> expectedUmiCounts = new HashMap<>();
		final TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(DigitalExpressionTest.EXPECTED_OUTFILE_LONG);
		for (final TabbedTextFileWithHeaderParser.Row row: parser) {
			final String cellBarcode = row.getField("CELL");
			final String geneName = row.getField("GENE");
			final int umiCount = row.getIntegerField("UMI_COUNT");
			expectedUmiCounts.put(new Pair<>(cellBarcode, geneName), umiCount);
		}
		Assert.assertEquals(umiCounts, expectedUmiCounts);
	}
}
