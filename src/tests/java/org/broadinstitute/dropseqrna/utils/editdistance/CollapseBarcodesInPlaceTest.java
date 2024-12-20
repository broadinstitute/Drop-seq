package org.broadinstitute.dropseqrna.utils.editdistance;

import java.io.File;
import java.io.IOException;
import java.util.*;

import org.broadinstitute.dropseqrna.utils.CompareBAMTagValues;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.testng.annotations.Test;

import org.testng.Assert;
import picard.nio.PicardHtsPath;

public class CollapseBarcodesInPlaceTest {


	private static File TEST_DATA = new File ("testdata/org/broadinstitute/dropseq/utils/N701_small.bam");
	// test substitutions with 100 reads core, 1 read non-core, read quality=10
	private static File TEST_SUB = new File ("testdata/org/broadinstitute/dropseq/utils/N701_small.expected_subsitutions.bam");
	// test indel+sub with 100 reads core, 1 read non-core, read quality=10
	private static File TEST_INDEL = new File ("testdata/org/broadinstitute/dropseq/utils/N701_small.expected_indels.bam");

	@Test
	public void filterBarcodesByNumReadsTest() {
		ObjectCounter<String> o = new ObjectCounter<>();
		o.incrementByCount("KEEP1", 5);
		o.incrementByCount("KEEP2", 4);
		o.incrementByCount("LOSE_AT_3", 2);
		o.incrementByCount("LOSE_AT_2", 1);

		CollapseBarcodesInPlace p = new CollapseBarcodesInPlace();
		ObjectCounter<String> actual = p.filterBarcodesByNumReads(o, 3);

		ObjectCounter<String> expected = new ObjectCounter<>();
		expected.incrementByCount("KEEP1", 5);
		expected.incrementByCount("KEEP2", 4);
		Assert.assertEquals(expected.getKeys(), actual.getKeys());

		actual = p.filterBarcodesByNumReads(o, 2);

		expected = new ObjectCounter<>();
		expected.incrementByCount("KEEP1", 5);
		expected.incrementByCount("KEEP2", 4);
		expected.incrementByCount("LOSE_AT_3", 2);
		Assert.assertEquals(expected.getKeys(), actual.getKeys());
	}

	@Test
	public void testCollapseBarcodes1() {
		ObjectCounter<String> o = new ObjectCounter<>();
		o.incrementByCount("TEST", 10);
		o.incrementByCount("FEST", 9); // TEST sub ed=1
		o.incrementByCount("ZEST", 8); // TEST sub ed=1
		o.incrementByCount("TEES", 7); //indel TEST ed=1
		o.incrementByCount("ZESZ", 6); //unrelated

		// test sub
		CollapseBarcodesInPlace p = new CollapseBarcodesInPlace();
		Map<String, String> actual = p.collapseBarcodes(o, false, 1);
		Map<String, String> expected = new HashMap<>();
		expected.put("FEST", "TEST");
		expected.put("ZEST", "TEST");
		Assert.assertEquals(expected, actual);

		// test alternate path 1 - all barcodes considered "core"
		List<String> barcodeList = o.getKeysOrderedByCount(true);
		actual = p.collapseBarcodes(barcodeList, o, false, 1);
		Assert.assertEquals(expected, actual);

		// test indel
		actual = p.collapseBarcodes(o, true, 1);
		expected = new HashMap<>();
		expected.put("FEST", "TEST");
		expected.put("ZEST", "TEST");
		expected.put("TEES", "TEST");
		Assert.assertEquals(expected, actual);

		// test alternate path - all barcodes considered "core"
		barcodeList = o.getKeysOrderedByCount(true);
		actual = p.collapseBarcodes(barcodeList, o, true, 1);
		Assert.assertEquals(expected, actual);
	}

	@Test
	public void testDoWorkSubstitution1 () {
		CollapseBarcodesInPlace p = new CollapseBarcodesInPlace();
		File out = getTempReportFile();
		out.deleteOnExit();
		List<File> files = new ArrayList<>();
		files.add(this.TEST_DATA);
		p.INPUT=files;
		p.PRIMARY_BARCODE="XC";
		p.FIND_INDELS=false;
		p.OUT_BARCODE="ZC";
		p.MIN_NUM_READS_CORE=100;
		// test data from unmapped BAM.
		p.MINIMUM_MAPPING_QUALITY=10;
		p.MIN_NUM_READS_NONCORE=1;
		// p.OUTPUT=new File ("testdata/org/broadinstitute/dropseq/utils/N701_small.expected_subsitutions.bam");
		// p.OUTPUT=new File ("testdata/org/broadinstitute/dropseq/utils/N701_small.expected_indels.bam");
		p.OUTPUT=out;
		p.doWork();

		List<String> tags = Collections.singletonList("ZC");

		// the two expected data sets differ.
		compareBAMTagValues(TEST_SUB, TEST_INDEL, tags, 1);

		// test against the correct answer for substitution
		compareBAMTagValues(TEST_SUB, out, tags, 0);

		// test against the wrong answer for substitution
		compareBAMTagValues(TEST_INDEL, out, tags, 1);
	}

	@Test
	public void testDoWorkSubstitution1Threaded () {
		CollapseBarcodesInPlace p = new CollapseBarcodesInPlace();
		p.NUM_THREADS=2;
		File out = getTempReportFile();
		out.deleteOnExit();
		List<File> files = new ArrayList<>();
		files.add(this.TEST_DATA);
		p.INPUT=files;
		p.PRIMARY_BARCODE="XC";
		p.FIND_INDELS=false;
		p.OUT_BARCODE="ZC";
		p.MIN_NUM_READS_CORE=100;
		// test data from unmapped BAM.
		p.MINIMUM_MAPPING_QUALITY=10;
		p.MIN_NUM_READS_NONCORE=1;
		p.OUTPUT=out;
		p.doWork();

		List<String> tags = Collections.singletonList("ZC");

		// the two expected data sets differ.
		compareBAMTagValues(TEST_SUB, TEST_INDEL, tags, 1);

		// test against the correct answer for substitution
		compareBAMTagValues(TEST_SUB, out, tags, 0);

		// test against the wrong answer for substitution
		compareBAMTagValues(TEST_INDEL, out, tags, 1);

	}
	private File getTempReportFile () {
		File tempFile=null;

		try {
			tempFile = File.createTempFile("CollapseBarcodesInPlaceTest", ".bam");
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		return tempFile;
	}

	private void compareBAMTagValues(File input1, File input2, List<String> tags, int expectedProgramValue) {
		CompareBAMTagValues cbtv = new CompareBAMTagValues();
		cbtv.INPUT_1 = Collections.singletonList(new PicardHtsPath(input1));
		cbtv.INPUT_2 = Collections.singletonList(new PicardHtsPath(input2));
		cbtv.TAGS_1 = tags;
		cbtv.TAGS_2 = tags;
		cbtv.STRICT = true;
		int result = cbtv.doWork();
		Assert.assertTrue(result == expectedProgramValue);
	}
}
