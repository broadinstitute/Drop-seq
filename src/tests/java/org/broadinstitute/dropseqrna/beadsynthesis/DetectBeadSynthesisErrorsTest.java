/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
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
