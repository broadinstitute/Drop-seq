package org.broadinstitute.dropseqrna.barnyard;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.StringWriter;
import java.util.Collections;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.nio.PicardHtsPath;

public class GatherMolecularBarcodeDistributionByGeneTest {

	private static final File IN_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene_retagged.bam");
	private static final File IN_CELL_BARCODE_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene.cellbarcodes.txt");
	private static final File OUT_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene_retagged.legacy_labels.molBC.txt");

	@Test
	// Note: This is a very basic test, because the output is basically the same as DGE testing, so we're covering the same ground.
	// The only difference is in the output of the read counts per UMI, and that's explicitly tested in the input vs output files here.
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
		g.INPUT = Collections.singletonList(new PicardHtsPath(IN_FILE));
		g.OUTPUT=outFile;
		g.LEGACY_COLUMN_LABELS = true;


        int result = g.doWork();
        Assert.assertEquals(result, 0);

        try {
			Assert.assertTrue (FileUtils.contentEquals(outFile, OUT_FILE));
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

	@Test
	public void testWritePerTranscriptStatsOrder() throws IOException {
		// Most supported UMI first, ties broken alphabetically, whatever order the UMIs were added in.
		final ObjectCounter<String> counts = new ObjectCounter<>();
		counts.incrementByCount("TTT", 2);
		counts.incrementByCount("CCC", 5);
		counts.incrementByCount("GGG", 2);
		counts.incrementByCount("AAA", 1);
		final StringWriter sw = new StringWriter();
		final BufferedWriter out = new BufferedWriter(sw);
		GatherMolecularBarcodeDistributionByGene.writePerTranscriptStats("G", "C", counts, out);
		out.flush();
		Assert.assertEquals(sw.toString().trim().split("\\R"), new String[]{"C\tG\tCCC\t5", "C\tG\tGGG\t2", "C\tG\tTTT\t2", "C\tG\tAAA\t1"});
	}
}
