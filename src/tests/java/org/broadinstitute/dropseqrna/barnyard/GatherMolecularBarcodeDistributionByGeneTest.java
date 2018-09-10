package org.broadinstitute.dropseqrna.barnyard;

import java.io.File;
import java.io.IOException;

import org.apache.commons.io.FileUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

public class GatherMolecularBarcodeDistributionByGeneTest {

	private static final File IN_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene_retagged.bam");
	private static final File IN_CELL_BARCODE_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene.cellbarcodes.txt");
	private static final File OUT_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene_retagged.molBC.txt");


	@Test
	public void testDoWork() {

		File outFile=null;
		try {
			outFile = File.createTempFile("testGatherMolecularBarcodeDistributionByGene.", ".digital_expression.txt");
	        outFile.deleteOnExit();
		} catch (IOException e) {
			e.printStackTrace();
		}

		GatherMolecularBarcodeDistributionByGene g = new GatherMolecularBarcodeDistributionByGene();
		g.CELL_BC_FILE=IN_CELL_BARCODE_FILE;
		g.INPUT=IN_FILE;
		g.OUTPUT=outFile;


        int result = g.doWork();
        Assert.assertEquals(result, 0);

        try {
			Assert.assertTrue (FileUtils.contentEquals(outFile, OUT_FILE));
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
