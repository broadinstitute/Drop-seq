package org.broadinstitute.dropseqrna.sbarro;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.io.FileUtils;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.broadinstitute.dropseqrna.utils.editdistance.MapBarcodesByEditDistance;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

public class BipartiteRabiesVirusCollapseTest {

	private static final File TEST_BAM = new File("testdata/org/broadinstitute/dropseq/sbarro/10_cells.bam");
	private static final File EXPECTED_REPORT_ED0 = new File("testdata/org/broadinstitute/dropseq/sbarro/10_cells.report.txt");
	
	// made the data smaller to speed up test time.
	private static final File TEST_BAM_SMALL = new File("testdata/org/broadinstitute/dropseq/sbarro/CTACCCAAGACCTAGG.bam");
	private static final File EXPECTED_REPORT_SMALL = new File("testdata/org/broadinstitute/dropseq/sbarro/CTACCCAAGACCTAGG.report.txt");
	
	
	@BeforeClass
	public void beforeSuite() {
		TestUtils.setInflaterDeflaterIfMacOs();
	}
	
	@Test
	public void testBipartiteUMISharing () throws IOException {
		BipartiteRabiesVirusCollapse c = new BipartiteRabiesVirusCollapse();
		c.OUT_TAG="rz";
		c.READ_MQ=0;
		c.MIN_COUNT=4;
		c.UMI_SHARING_ONLY_MODE=true;
		c.UMI_SHARING_EDIT_DISTANCE=0;
		c.UMI_SHARING_THRESHOLD=0.8;
		c.INPUT=TEST_BAM;
		c.OUTPUT=File.createTempFile("BipartiteRabiesVirusCollapseTest", ".bam");
		c.OUTPUT.deleteOnExit();
		c.REPORT=File.createTempFile("BipartiteRabiesVirusCollapseTest", ".report.txt");
		c.REPORT.deleteOnExit();
		c.doWork();
		
		Assert.assertTrue (FileUtils.contentEquals(c.REPORT, EXPECTED_REPORT_ED0));
		
		
	}
	
	@Test (enabled = true)
	public void testBipartiteWithoutUMISharing () throws IOException {
		BipartiteRabiesVirusCollapse c = new BipartiteRabiesVirusCollapse();
		c.OUT_TAG="rz";
		c.READ_MQ=0;
		c.MIN_COUNT=4;
		c.UMI_SHARING_ONLY_MODE=false;
		c.INPUT=TEST_BAM_SMALL;
		c.OUTPUT=File.createTempFile("BipartiteRabiesVirusCollapseTest", ".bam");
		c.OUTPUT.deleteOnExit();
		c.REPORT=File.createTempFile("BipartiteRabiesVirusCollapseTest", ".no_sharing.report.txt");
		c.REPORT.deleteOnExit();
		c.doWork();
		
		Assert.assertTrue (FileUtils.contentEquals(c.REPORT, EXPECTED_REPORT_SMALL));
		
		
	}
	
	
	@Test
	public void collapseRabiesBarcodesEx1() {

		BipartiteRabiesVirusCollapse c = new BipartiteRabiesVirusCollapse();
		MapBarcodesByEditDistance med = new MapBarcodesByEditDistance(false);

		ObjectCounter<String> t1 = new ObjectCounter<>();
		t1.incrementByCount("TCATAATCGCGTCTGTAATT", 1117);
		t1.incrementByCount("CCCCTAGCAGGTCCGTGGGA", 961);
		t1.incrementByCount("TACTAAAACCGTCCGTGGGA", 4);
		t1.incrementByCount("CCCCTANNAGGTCCGNNNNN", 3);

		Map<String, BipartiteRabiesVirusCollapseResultCollection> result = c.collapseRabiesBarcodes(t1, 10, false, 1,
				med);
		Map<String, Collection<String>> actual = result.get("TACTAAAACCGTCCGTGGGA").getBarcodeMapping();
		Map<String, Collection<String>> expected = new HashMap<>();
		expected.put("TACTAAAACCGTCCGTGGGA", Arrays.asList("CCCCTAGCAGGTCCGTGGGA"));
		Assert.assertEquals(actual, expected);
	}

	@Test
	public void collapseRabiesBarcodesEx2() {

		BipartiteRabiesVirusCollapse c = new BipartiteRabiesVirusCollapse();
		MapBarcodesByEditDistance med = new MapBarcodesByEditDistance(false);

		ObjectCounter<String> t1 = new ObjectCounter<>();
		t1.incrementByCount("AAGAACACTTGAGAAGTGCT", 100);
		t1.incrementByCount("TCTAACGTATTCCTGTACGG", 50);
		t1.incrementByCount("AAGAACACTTTCCTGTACGG", 1);

		Map<String, BipartiteRabiesVirusCollapseResultCollection> result = c.collapseRabiesBarcodes(t1, 10, false, 1,
				med);
		String f = result.get("AAGAACACTTTCCTGTACGG").toString();
		Map<String, Collection<String>> actual = result.get("AAGAACACTTTCCTGTACGG").getBarcodeMapping();
		Map<String, Collection<String>> expected = new HashMap<>();
		expected.put("AAGAACACTTTCCTGTACGG", Arrays.asList("AAGAACACTTGAGAAGTGCT", "TCTAACGTATTCCTGTACGG"));
		Assert.assertEquals(actual, expected);
	}

	@Test
	public void collapseRabiesBarcodesEx3() {

		BipartiteRabiesVirusCollapse c = new BipartiteRabiesVirusCollapse();
		MapBarcodesByEditDistance med = new MapBarcodesByEditDistance(false);

		ObjectCounter<String> t1 = new ObjectCounter<>();
		t1.incrementByCount("TCGATCCGACAACTATAGCG", 100);
		t1.incrementByCount("CCCTCCCTTGCGCGCATCGC", 50);
		t1.incrementByCount("AAGTAGCCAACCTGGGTGGA", 1);
		t1.incrementByCount("TCTATAATACAACTATAGCG", 1);
		t1.incrementByCount("CAGCACGCATTGTTAAATCT", 1);

		Map<String, BipartiteRabiesVirusCollapseResultCollection> result = c.collapseRabiesBarcodes(t1, 10, false, 1,
				med);
		Map<String, Collection<String>> actual = result.get("TCTATAATACAACTATAGCG").getBarcodeMapping();
		Map<String, Collection<String>> expected = new HashMap<>();
		expected.put("TCTATAATACAACTATAGCG", Arrays.asList("TCGATCCGACAACTATAGCG"));
		Assert.assertEquals(actual, expected);
	}

}
