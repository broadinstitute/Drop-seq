package org.broadinstitute.dropseqrna.spermseq.metrics.spermalleles;

import java.io.File;
import java.io.IOException;

import org.apache.commons.io.FileUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

public class GenotypeSpermTest {
	
	private static final File BAM_FILE = new File("testdata/org/broadinstitute/spermseq/spermalleles/GenotypeSperm.bam");
	private static final File CELL_BARCODE_FILE = new File("testdata/org/broadinstitute/spermseq/spermalleles/GenotypeSperm.cellBarcodes.txt");
	private static final File INTERVALS_FILE = new File("testdata/org/broadinstitute/spermseq/spermalleles/GenotypeSperm.intervals");
	private static final File EXPECTED_RESULTS_FILE = new File("testdata/org/broadinstitute/spermseq/spermalleles/GenotypeSperm.result.txt");
	
  @Test
  public void testFullProgram() {
	  	File outFile=null;
		try {
			outFile = File.createTempFile("GenotypeSpermTest.", ".result.txt");
	        outFile.deleteOnExit();
		} catch (IOException e) {
			e.printStackTrace();
		}
		GenotypeSperm dsa = new GenotypeSperm();
		dsa.INPUT=BAM_FILE;
		dsa.CELL_BC_FILE=CELL_BARCODE_FILE;
		dsa.INTERVALS=INTERVALS_FILE;
		dsa.OUTPUT=outFile;
		      
      int result = dsa.doWork();
      Assert.assertEquals(result, 0);
      try {
			Assert.assertTrue (FileUtils.contentEquals(outFile, EXPECTED_RESULTS_FILE));
		} catch (IOException e) {
			e.printStackTrace();
		}
  }
  
}
