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

import org.broadinstitute.dropseqrna.utils.BaseRange;
import org.testng.Assert;
import org.testng.annotations.Test;

public class FilterValidRabiesBarcodesTest {


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
