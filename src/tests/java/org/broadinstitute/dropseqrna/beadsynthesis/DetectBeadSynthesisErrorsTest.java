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

import org.junit.Assert;
import org.testng.annotations.Test;


public class DetectBeadSynthesisErrorsTest {

	DetectBeadSynthesisErrors gbse = new DetectBeadSynthesisErrors();
	private final String cellBCTag = gbse.CELL_BARCODE_TAG;
	private final String molBCTag = gbse.MOLECULAR_BARCODE_TAG;
	
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
	
}
