package org.broadinstitute.dropseqrna.beadsynthesis;

import java.util.Map;

import htsjdk.samtools.SAMRecord;

import org.junit.Assert;
import org.testng.annotations.Test;


public class GatherBeadSynthesisErrorsTest {

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
	
	public void padCellBarcodeTest2() {
		DetectBeadSynthesisErrors gbse = new DetectBeadSynthesisErrors();
		String startBarcode = "AAAAAAAAAAAZ";
		int errorBase = 8;
		int umiLength = 8;
		
		String fixedBarcode=gbse.padCellBarcode(startBarcode, errorBase, umiLength);
		String expected="AAAAAAAAAAAN";
		Assert.assertEquals(expected, fixedBarcode);
		
		
	}
	
	public void padCellBarcodeTest3() {
		DetectBeadSynthesisErrors gbse = new DetectBeadSynthesisErrors();
		String startBarcode = "AAAAAAAAAAAZ";
		int errorBase = -1;
		int umiLength = 8;
		
		String fixedBarcode=gbse.padCellBarcode(startBarcode, errorBase, umiLength);
		String expected="AAAAAAAAAAAZ";
		Assert.assertEquals(expected, fixedBarcode);		
	}
	
}
