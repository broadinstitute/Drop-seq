package org.broadinstitute.dropseqrna.barnyard;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.dropseqrna.beadsynthesis.GenerateRandomUMIs;
import org.testng.annotations.Test;

import junit.framework.Assert;

public class SelectCellsByNumTranscriptsTest {

	File SINGLE_ORGANISM_BAM= new File ("testdata/org/broadinstitute/dropseq/utils/N701_small.bam");
	File SINGLE_ORGANISM_BAM_EXPECTED_CB_100_READS=new File ("testdata/org/broadinstitute/dropseq/utils/N701_small.cell_barcodes_100_reads.txt");
	File SINGLE_ORGANISM_BAM_EXPECTED_CB_100_TRANSCRIPTS=new File ("testdata/org/broadinstitute/dropseq/utils/N701_small.cell_barcodes_100_transcripts.txt");

	@Test
	public void testDoWorkSingleOrganism() {
		SelectCellsByNumTranscripts s = new SelectCellsByNumTranscripts();
		File outFile = getTempReportFile("SelectCellsByNumTranscripts.", ".cellBarcodes");
		File outFileMinReads = getTempReportFile("SelectCellsByNumTranscripts.", ".minReads");
		outFile.deleteOnExit();
		outFileMinReads.deleteOnExit();

		s.INPUT=this.SINGLE_ORGANISM_BAM;
		s.MIN_TRANSCRIPTS_PER_CELL=100;
		s.OUTPUT=outFile;
		s.MIN_READS_PER_CELL=100;
		s.OUTPUT_INTERIM_CELLS=outFileMinReads;
		s.READ_MQ=0;
		s.ORGANISM=Collections.emptyList();
		int success = s.doWork();

		try {
			boolean t1 = FileUtils.contentEquals(outFileMinReads, SINGLE_ORGANISM_BAM_EXPECTED_CB_100_READS);
			boolean t2 = FileUtils.contentEquals(outFile, SINGLE_ORGANISM_BAM_EXPECTED_CB_100_TRANSCRIPTS);
			Assert.assertTrue(t1);
			Assert.assertTrue(t2);
		} catch (IOException e) {
			e.printStackTrace();
		}

		Assert.assertTrue(success==0);
	}

	@Test
	public void testReadWriteBarcodes () {
		GenerateRandomUMIs g = new GenerateRandomUMIs(0);
		List<String> barcodes = new ArrayList<>();
		for (int i=0; i<20; i++)
			barcodes.add(g.getRandomString(12));
		File outFile = getTempReportFile("SelectCellsByNumTranscripts.", ".cellBarcodes");
		outFile.deleteOnExit();
		SelectCellsByNumTranscripts.writeBarcodes(outFile, barcodes);
		List<String> barcodes2 = SelectCellsByNumTranscripts.readBarcodes(outFile);
		Assert.assertEquals(barcodes, barcodes2);

	}

	// Moar code coverage.
	@Test (expectedExceptions = htsjdk.samtools.SAMException.class)
	public void testFileWriteException () {
		SelectCellsByNumTranscripts.writeBarcodes(new File("/foo/bar/zoo"), null);
	}

	@Test (expectedExceptions = htsjdk.samtools.SAMException.class)
	public void testFileReadException () {
		SelectCellsByNumTranscripts.readBarcodes(new File("/foo/bar/zoo"));
	}

	private File getTempReportFile (final String prefix, final String suffix) {
		File tempFile=null;

		try {
			tempFile = File.createTempFile(prefix, suffix);
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		return tempFile;
	}
}
