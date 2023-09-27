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
package org.broadinstitute.dropseqrna.beadsynthesis;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.dropseqrna.TranscriptomeException;
import org.broadinstitute.dropseqrna.utils.CompareBAMTagValues;
import org.testng.Assert;
import org.testng.annotations.Test;


public class DetectBeadSynthesisErrorsTest {

	// test indel+sub with 100 reads core, 1 read non-core, read quality=10
	private static File TEST_INDEL = new File ("testdata/org/broadinstitute/dropseq/utils/N701_small.expected_indels.bam");
	private static File EXPECTED_REPORT = new File ("testdata/org/broadinstitute/dropseq/beadsynthesis/DetectBeadSynthesisErrors.report");
	private static File EXPECTED_STATS = new File ("testdata/org/broadinstitute/dropseq/beadsynthesis/DetectBeadSynthesisErrors.stats");
	private static File EXPECTED_SUMMARY = new File ("testdata/org/broadinstitute/dropseq/beadsynthesis/DetectBeadSynthesisErrors.summary");
	private static File EXPECTED_BAM = new File ("testdata/org/broadinstitute/dropseq/beadsynthesis/DetectBeadSynthesisErrors.bam");

		@Test
	public void testDoWork() {
		DetectBeadSynthesisErrors gbse = new DetectBeadSynthesisErrors();

		File report = getTempReportFile("DetectBeadSynthesisErrorsTest", ".report");
		File stats = getTempReportFile("DetectBeadSynthesisErrorsTest", ".stats");
		File summary = getTempReportFile("DetectBeadSynthesisErrorsTest", ".summary");
		File cleanBAM = getTempReportFile("DetectBeadSynthesisErrorsTest", ".bam");

		report.deleteOnExit();
		stats.deleteOnExit();
		summary.deleteOnExit();
		cleanBAM.deleteOnExit();

		gbse.INPUT=Arrays.asList(TEST_INDEL);
		gbse.NUM_BARCODES=1000; // more than the BAM has.
		gbse.SUMMARY=summary;
	    gbse.OUTPUT=cleanBAM;
		gbse.REPORT=report;
		gbse.OUTPUT_STATS=stats;

		// test custom command line validation
		gbse.MAX_BARCODE_ERRORS_IN_RAM=5;  // i want to test serialization/deserialization.

		int result = gbse.doWork();
		Assert.assertEquals(0, result);

		try {
			boolean t1 = FileUtils.contentEquals(report, EXPECTED_REPORT);
			boolean t2 = FileUtils.contentEquals(stats, EXPECTED_STATS);
			boolean t3 = FileUtils.contentEquals(summary, EXPECTED_SUMMARY);
			Assert.assertTrue(t1);
			Assert.assertTrue(t2);
			Assert.assertTrue(t3);
		} catch (IOException e) {
			e.printStackTrace();
		}

		// test BAM
		CompareBAMTagValues cbtv = new CompareBAMTagValues();
		cbtv.INPUT_1=EXPECTED_BAM;
		cbtv.INPUT_2=cleanBAM;
		List<String> tags = new ArrayList<>();
		tags.add("XC");
		cbtv.TAGS=tags;
		int r = cbtv.doWork();
		Assert.assertTrue(r==0);




	}


	@Test
	public void padCellBarcodeTest1() {
		DetectBeadSynthesisErrors gbse = new DetectBeadSynthesisErrors();
		String startBarcode = "AAAAAAAAAAZZ";
		int errorBase = 7;
		int umiLength = 8;

		String fixedBarcode=gbse.padCellBarcode(startBarcode, errorBase, umiLength);
		String expected="AAAAAAAAAANN";
		Assert.assertEquals(expected, fixedBarcode);
	}

	@Test
	public void padCellBarcodeTest2() {
		DetectBeadSynthesisErrors gbse = new DetectBeadSynthesisErrors();
		String startBarcode = "AAAAAAAAAAAZ";
		int errorBase = 8;
		int umiLength = 8;

		String fixedBarcode=gbse.padCellBarcode(startBarcode, errorBase, umiLength);
		String expected="AAAAAAAAAAAN";
		Assert.assertEquals(expected, fixedBarcode);


	}
	@Test
	public void padCellBarcodeTest3() {
		DetectBeadSynthesisErrors gbse = new DetectBeadSynthesisErrors();
		String startBarcode = "AAAAAAAAAAAZ";
		int errorBase = -1;
		int umiLength = 8;

		String fixedBarcode=gbse.padCellBarcode(startBarcode, errorBase, umiLength);
		String expected="AAAAAAAAAAAZ";
		Assert.assertEquals(expected, fixedBarcode);
	}

	@Test
	public void fixUMITest1 () {
		DetectBeadSynthesisErrors gbse = new DetectBeadSynthesisErrors();
		String startBarcode = "AAAAAAAAAAAZ";
		String umi="GGGGGGGT";
		int errorBase = 8;
		int umiLength = 8;

		String fixedUMI=gbse.fixUMI(startBarcode, umi, errorBase);
		String expected="ZGGGGGGG";
		Assert.assertEquals(expected, fixedUMI);
	}

	@Test
	public void fixUMITest2 () {
		DetectBeadSynthesisErrors gbse = new DetectBeadSynthesisErrors();
		String startBarcode = "ACGCTCATACAG";
		String umi="TCCTTATT";
		int errorBase = 7;

		String fixedUMI=gbse.fixUMI(startBarcode, umi, errorBase);
		String expected="AGTCCTTA";
		Assert.assertEquals(expected, fixedUMI);

		String fixedBarcode=gbse.padCellBarcode(startBarcode, errorBase, umi.length());
		String expectedCell="ACGCTCATACNN";
		Assert.assertEquals(expectedCell, fixedBarcode);

	}

	@Test
	public void fixUMITest3 () {
		DetectBeadSynthesisErrors gbse = new DetectBeadSynthesisErrors();
		String startBarcode = "GAGCTAGTTACT";
		String umi="ATCTTTTT";
		int errorBase = 7;

		String fixedUMI=gbse.fixUMI(startBarcode, umi, errorBase);
		String expected="CTATCTTT";
		Assert.assertEquals(expected, fixedUMI);

		String fixedBarcode=gbse.padCellBarcode(startBarcode, errorBase, umi.length());
		String expectedCell="GAGCTAGTTANN";
		Assert.assertEquals(expectedCell, fixedBarcode);

	}


	@Test
	public void testBuildBarcodeNeighborGroups () {
		List<BeadSynthesisErrorData> d= new ArrayList<>();
		double umiBiasThreshold=0.8;

		// add a set of 4 neighbors.
		d.add(generateBaseCounts("AAAA1", BeadSynthesisErrorType.SYNTH_MISSING_BASE, 100, 9, umiBiasThreshold));
		d.add(generateBaseCounts("AAAA2", BeadSynthesisErrorType.SYNTH_MISSING_BASE, 100, 9, umiBiasThreshold));
		d.add(generateBaseCounts("AAAA3", BeadSynthesisErrorType.SYNTH_MISSING_BASE, 100, 9, umiBiasThreshold));
		d.add(generateBaseCounts("AAAA4", BeadSynthesisErrorType.SYNTH_MISSING_BASE, 100, 9, umiBiasThreshold));

		// add a set of 3 neighbors
		d.add(generateBaseCounts("BBBB1", BeadSynthesisErrorType.SYNTH_MISSING_BASE, 100, 9, umiBiasThreshold));
		d.add(generateBaseCounts("BBBB2", BeadSynthesisErrorType.SYNTH_MISSING_BASE, 100, 9, umiBiasThreshold));
		d.add(generateBaseCounts("BBBB3", BeadSynthesisErrorType.SYNTH_MISSING_BASE, 100, 9, umiBiasThreshold));

		// add a set of 2 neighbors
		d.add(generateBaseCounts("CCCC1", BeadSynthesisErrorType.SYNTH_MISSING_BASE, 100, 9, umiBiasThreshold));
		d.add(generateBaseCounts("CCCC2", BeadSynthesisErrorType.SYNTH_MISSING_BASE, 100, 9, umiBiasThreshold));

		// add a set of 1 neighbor
		d.add(generateBaseCounts("DDDD1", BeadSynthesisErrorType.SYNTH_MISSING_BASE, 100, 9, umiBiasThreshold));

		// add a set of 2 neighbors with errors and 2 without errors. Only the 2 errors should be returned.
		d.add(generateBaseCounts("EEEE1", BeadSynthesisErrorType.SYNTH_MISSING_BASE, 100, 9, umiBiasThreshold));
		d.add(generateBaseCounts("EEEE2", BeadSynthesisErrorType.SYNTH_MISSING_BASE, 100, 9, umiBiasThreshold));
		d.add(generateBaseCounts("EEEE3", BeadSynthesisErrorType.NO_ERROR, 100, 9, umiBiasThreshold));
		d.add(generateBaseCounts("EEEE4", BeadSynthesisErrorType.NO_ERROR, 100, 9, umiBiasThreshold));

		DetectBeadSynthesisErrors dbse = new DetectBeadSynthesisErrors();
		Collection <BarcodeNeighborGroup> result = dbse.buildBarcodeNeighborGroups(d, umiBiasThreshold).values();

		for (BarcodeNeighborGroup g: result) {
			if (g.getRootSequence().equals("AAAA"))
				Assert.assertEquals(4, g.getNeighborCellBarcodes().size());
			if (g.getRootSequence().equals("BBBB"))
				Assert.assertEquals(3, g.getNeighborCellBarcodes().size());
			if (g.getRootSequence().equals("CCCC"))
				Assert.assertEquals(2, g.getNeighborCellBarcodes().size());
			if (g.getRootSequence().equals("DDDD"))
				Assert.assertEquals(1, g.getNeighborCellBarcodes().size());
			if (g.getRootSequence().equals("EEEE"))
				Assert.assertEquals(2, g.getNeighborCellBarcodes().size());
		}
	}

	// methods for generating BaseDistributionMetricCollection with a few patterns.
	private BeadSynthesisErrorData generateBaseCounts (final String cellBarcode, final BeadSynthesisErrorType error, final int numUMis, final int numUMIBases, final double umiBiasThreshold) {
		GenerateRandomUMIs gru = new GenerateRandomUMIs(umiBiasThreshold);
		Collection<String> umis = gru.getUMICollection(numUMis, numUMIBases, error);
		BeadSynthesisErrorData d = new BeadSynthesisErrorData(cellBarcode);
		d.addUMI(umis);
		// validate that the generated data has the right error type.
		BeadSynthesisErrorType actualError = d.getErrorType(umiBiasThreshold);
		Assert.assertEquals(error, actualError);
		return d;
	}

	private File getTempReportFile (final String prefix, final String suffix) {
		File tempFile=null;

		try {
			tempFile = File.createTempFile(prefix, suffix);
			tempFile.deleteOnExit();
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		return tempFile;
	}

	@Test
	public void testCustomCommandParsing () {
		DetectBeadSynthesisErrors gbse = new DetectBeadSynthesisErrors();
		gbse.INPUT=Arrays.asList(TEST_INDEL);
		gbse.OUTPUT_STATS=getTempReportFile("DetectBeadSynthesisErrorsTest", ".stats");
		gbse.SUMMARY=getTempReportFile("DetectBeadSynthesisErrorsTest", ".summary");

		gbse.CELL_BC_FILE=null;
		gbse.NUM_BARCODES=2;
		String [] errors = gbse.customCommandLineValidation();
		Assert.assertNull(errors);

		gbse.CELL_BC_FILE=EXPECTED_REPORT;
		gbse.NUM_BARCODES=null;
		errors = gbse.customCommandLineValidation();
		Assert.assertNull(errors);

	}

	@Test(expectedExceptions=TranscriptomeException.class)
	public void testCustomCommandParsing2 () {
		DetectBeadSynthesisErrors gbse = new DetectBeadSynthesisErrors();
		gbse.INPUT=Arrays.asList(TEST_INDEL);
		gbse.OUTPUT_STATS=getTempReportFile("DetectBeadSynthesisErrorsTest", ".stats");
		gbse.SUMMARY=getTempReportFile("DetectBeadSynthesisErrorsTest", ".summary");

		gbse.CELL_BC_FILE=EXPECTED_REPORT;
		gbse.NUM_BARCODES=2;
		String [] errors = gbse.customCommandLineValidation();
		Assert.assertTrue(errors.length>0);
	}



}
