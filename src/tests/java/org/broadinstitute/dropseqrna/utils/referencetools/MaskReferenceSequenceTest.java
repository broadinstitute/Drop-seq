package org.broadinstitute.dropseqrna.utils.referencetools;

import java.io.File;
import java.io.IOException;

import org.apache.commons.io.FileUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

public class MaskReferenceSequenceTest {

	private static final File IN_REF = new File("testdata/org/broadinstitute/dropseq/utils/referencetools/fake_ref.fasta");
	private static final File OUT_REF_INTERVALS = new File("testdata/org/broadinstitute/dropseq/utils/referencetools/fake_ref.filtered_by_intervals.fasta");
	private static final File OUT_REF_CONTIGS = new File("testdata/org/broadinstitute/dropseq/utils/referencetools/fake_ref.filtered_by_contigs.fasta");
	private static final File INTERVAL_FILE = new File("testdata/org/broadinstitute/dropseq/utils/referencetools/fake_ref.intervals");

	@Test
	public void testDoWork() {
		File outFile=null;
		try {
			outFile = File.createTempFile("MaskReferenceSequenceTest.", ".fasta");
			outFile.deleteOnExit();
		} catch (IOException e) {
			e.printStackTrace();
		}

		MaskReferenceSequence m = new MaskReferenceSequence();
		String [] args = new String [4];
		args[0]="OUTPUT="+outFile.getAbsolutePath();
		args[1]="INTERVAL_LIST="+INTERVAL_FILE.getAbsolutePath();
		args[2]="REFERENCE_SEQUENCE="+IN_REF.getAbsolutePath();
		args[3]="OUTPUT_LINE_LENGTH=50";

		int r = m.instanceMain(args);
		Assert.assertTrue(r==0);

		try {
			Assert.assertTrue (FileUtils.contentEquals(outFile, OUT_REF_INTERVALS));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	@Test
	public void testDoWork2() {
		File outFile=null;
		try {
			outFile = File.createTempFile("MaskReferenceSequenceTest.", ".fasta");
			outFile.deleteOnExit();
		} catch (IOException e) {
			e.printStackTrace();
		}

		MaskReferenceSequence m = new MaskReferenceSequence();
		String [] args = new String [4];
		args[0]="OUTPUT="+outFile.getAbsolutePath();
		args[1]="CONTIG_PATTERN_TO_IGNORE=fake_contig_2";
		args[2]="REFERENCE_SEQUENCE="+IN_REF.getAbsolutePath();
		args[3]="OUTPUT_LINE_LENGTH=50";
		int r = m.instanceMain(args);
		Assert.assertTrue(r==0);
		try {
			Assert.assertTrue (FileUtils.contentEquals(outFile, OUT_REF_CONTIGS));
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
