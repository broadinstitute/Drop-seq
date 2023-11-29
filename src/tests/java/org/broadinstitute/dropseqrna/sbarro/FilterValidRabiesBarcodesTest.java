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
package org.broadinstitute.dropseqrna.sbarro;

import java.io.File;
import java.io.IOException;

import org.broadinstitute.dropseqrna.utils.BaseRange;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

public class FilterValidRabiesBarcodesTest {


	private static final File INPUT = new File("testdata/org/broadinstitute/dropseq/sbarro/TagReadWithRabiesBarcodes_output.bam");
	private static final File EXPECTED_OUTPUT_ACCEPTED = new File("testdata/org/broadinstitute/dropseq/sbarro/FilterValidRabiesBarcodes.output.bam");
	private static final File EXPECTED_OUTPUT_REJECTED = new File("testdata/org/broadinstitute/dropseq/sbarro/FilterValidRabiesBarcodes.output_rejectedd.bam");
	
	
	@Test
	public void endToEnd () throws IOException {
		FilterValidRabiesBarcodes f = new FilterValidRabiesBarcodes();
		
		f.INPUT=INPUT;
		f.OUTPUT_ACCEPTED=File.createTempFile("FilterValidRabiesBarcodesTest", ".accepted.bam");
		f.OUTPUT_REJECTED=File.createTempFile("FilterValidRabiesBarcodesTest", ".rejected.bam");
		f.OUTPUT_ACCEPTED.deleteOnExit();
		f.OUTPUT_REJECTED.deleteOnExit();
		f.STOP_BARCODE_LENGTH_RANGE="10-10";
		f.POLYA_BARCODE_LENGTH_RANGE="10-10";
		f.MAX_GFP_ANCHOR_EDIT_DISTANCE=5;
		f.MAX_CASSETTE_ANCHOR_EDIT_DISTANCE=5;
		f.MAX_N_IN_BARCODE=0;
		
		int r = f.doWork();
		Assert.assertEquals(r, 0);
		
		TestUtils.assertSamFilesSame(EXPECTED_OUTPUT_ACCEPTED, f.OUTPUT_ACCEPTED, false);
		TestUtils.assertSamFilesSame(EXPECTED_OUTPUT_REJECTED, f.OUTPUT_REJECTED, false);
		
		
	}
	
	@Test
	// test with all fields not-null.
	public void testSimplePass () {

		Integer readGFPAnchorED=1;
		Integer readCassetteAnchorED=0;
		Integer stopCodonBCSize=11;
		Integer polyABCSize=12;

		BaseRange stopBCRange = new BaseRange(10, 12);
		BaseRange polyABCRange = new BaseRange(10, 12);
		Integer maxGFPED=1;
		Integer maxCassetteED=1;

		boolean result = FilterValidRabiesBarcodes.acceptRead(readGFPAnchorED, readCassetteAnchorED, stopCodonBCSize, polyABCSize, stopBCRange, polyABCRange, maxGFPED, maxCassetteED);
		Assert.assertTrue(result);

		// iteratively pass in null values for filters.  Removing filters doesn't change the result.
		result = FilterValidRabiesBarcodes.acceptRead(readGFPAnchorED, readCassetteAnchorED, stopCodonBCSize, polyABCSize, null, polyABCRange, maxGFPED, maxCassetteED);
		Assert.assertTrue(result);

		result = FilterValidRabiesBarcodes.acceptRead(readGFPAnchorED, readCassetteAnchorED, stopCodonBCSize, polyABCSize, stopBCRange, null, maxGFPED, maxCassetteED);
		Assert.assertTrue(result);

		result = FilterValidRabiesBarcodes.acceptRead(readGFPAnchorED, readCassetteAnchorED, stopCodonBCSize, polyABCSize, stopBCRange, polyABCRange, null, maxCassetteED);
		Assert.assertTrue(result);

		result = FilterValidRabiesBarcodes.acceptRead(readGFPAnchorED, readCassetteAnchorED, stopCodonBCSize, polyABCSize, stopBCRange, polyABCRange, maxGFPED, null);
		Assert.assertTrue(result);
	}

	@Test
	public void testSimpleReject() {
		Integer readGFPAnchorED=2;
		Integer readCassetteAnchorED=0;
		Integer stopCodonBCSize=11;
		Integer polyABCSize=10;

		BaseRange stopBCRange = new BaseRange(10, 10);
		BaseRange polyABCRange = new BaseRange(10, 10);
		Integer maxGFPED=0;
		Integer maxCassetteED=0;

		// initial test rejects because GFP anchor ED too big, stopCodonBC size too big.
		boolean result = FilterValidRabiesBarcodes.acceptRead(readGFPAnchorED, readCassetteAnchorED, stopCodonBCSize, polyABCSize, stopBCRange, polyABCRange, maxGFPED, maxCassetteED);
		Assert.assertFalse(result);

		// remove filters for GFP anchor ED and stopCodonBC size, should pass
		result = FilterValidRabiesBarcodes.acceptRead(readGFPAnchorED, readCassetteAnchorED, stopCodonBCSize, polyABCSize, null, polyABCRange, null, maxCassetteED);
		Assert.assertTrue(result);

		// remove just one filter for GFP anchor ED and stopCodonBC size, should fail
		result = FilterValidRabiesBarcodes.acceptRead(readGFPAnchorED, readCassetteAnchorED, stopCodonBCSize, polyABCSize, stopBCRange, polyABCRange, null, maxCassetteED);
		Assert.assertFalse(result);

		// remove just one filter for GFP anchor ED and stopCodonBC size, should fail
		result = FilterValidRabiesBarcodes.acceptRead(readGFPAnchorED, readCassetteAnchorED, stopCodonBCSize, polyABCSize, null, polyABCRange, maxGFPED, maxCassetteED);
		Assert.assertFalse(result);
	}
	
	
}
