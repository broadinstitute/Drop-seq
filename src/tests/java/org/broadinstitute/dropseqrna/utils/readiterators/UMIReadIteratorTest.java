package org.broadinstitute.dropseqrna.utils.readiterators;

import java.io.File;
import java.util.Collections;
import java.util.List;

import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMRecord;

public class UMIReadIteratorTest {

	private static final File INPUT = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene_retagged.bam");
	private static final File IN_CELL_BARCODE_FILE = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene.cellbarcodes.txt");

	@Test
	public void testIterate() {
		SamHeaderAndIterator headerAndIterator = SamFileMergeUtil.mergeInputs(Collections.singletonList(this.INPUT), false);
		List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(IN_CELL_BARCODE_FILE);

		UMIReadIterator i = new UMIReadIterator(headerAndIterator, "GE", "XC", "XM", "GS", true, 10, true, cellBarcodes);
		while (i.hasNext()) {
			List<SAMRecord> r = i.next();
			Assert.assertNotNull(r);
			Assert.assertTrue(r.size()>0);
		}

	}
}
