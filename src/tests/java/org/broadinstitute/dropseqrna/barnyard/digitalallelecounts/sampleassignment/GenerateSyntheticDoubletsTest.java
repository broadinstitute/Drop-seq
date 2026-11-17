package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.testng.Assert;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class GenerateSyntheticDoubletsTest {

	private static final File INPUT = new File(
			"testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts/sampleassignment/small_data_retagged.sam");
	// The input has 4 reads for each of these two cells.
	private static final String CELL_A = "ATCAGGGACAGA";
	private static final String CELL_B = "TTGCCTTACGCG";

	private GenerateSyntheticDoublets makeBase() {
		GenerateSyntheticDoublets g = new GenerateSyntheticDoublets();
		g.INPUT = Collections.singletonList(INPUT);
		g.OUTPUT = TestUtils.getTempReportFile("GenerateSyntheticDoublets.", ".bam");
		g.SUMMARY = TestUtils.getTempReportFile("GenerateSyntheticDoublets.", ".summary.txt");
		return g;
	}

	private static Map<String, Integer> countByTag(final File bam) throws IOException {
		Map<String, Integer> counts = new TreeMap<>();
		try (SamReader reader = SamReaderFactory.makeDefault().open(bam)) {
			for (SAMRecord r : reader)
				counts.merge(r.getStringAttribute("XC"), 1, Integer::sum);
		}
		return counts;
	}

	@Test
	public void testMultiCellFile() throws IOException {
		final File multiCellFile = File.createTempFile("multi_cell", ".txt");
		multiCellFile.deleteOnExit();
		Files.write(multiCellFile.toPath(), Collections.singletonList(CELL_A + ":" + CELL_B));
		GenerateSyntheticDoublets g = makeBase();
		g.MULTI_CELL_FILE = multiCellFile;
		Assert.assertEquals(g.doWork(), 0);

		// Only the merged cell is emitted, and it holds the reads of both cells.
		final List<String> summary = Files.readAllLines(g.SUMMARY.toPath());
		Assert.assertEquals(summary.size(), 1);
		Map<String, Integer> counts = countByTag(g.OUTPUT);
		Assert.assertEquals(counts.size(), 1);
		Assert.assertEquals(counts.get(summary.get(0)), Integer.valueOf(8));
		Assert.assertTrue(summary.get(0).contains(CELL_A) && summary.get(0).contains(CELL_B));
	}

	@Test
	public void testEmitSingletons() throws IOException {
		final File multiCellFile = File.createTempFile("multi_cell", ".txt");
		multiCellFile.deleteOnExit();
		Files.write(multiCellFile.toPath(), Collections.singletonList(CELL_A + ":" + CELL_B));
		GenerateSyntheticDoublets g = makeBase();
		g.MULTI_CELL_FILE = multiCellFile;
		g.EMIT_SINGLETONS = true;
		Assert.assertEquals(g.doWork(), 0);

		Map<String, Integer> counts = countByTag(g.OUTPUT);
		// Original reads for both cells plus 8 reads for the merged cell.
		Assert.assertEquals(counts.size(), 3);
		Assert.assertEquals(counts.get(CELL_A), Integer.valueOf(4));
		Assert.assertEquals(counts.get(CELL_B), Integer.valueOf(4));
	}

	@Test
	public void testRandomMultiCellsFromBarcodeFile() throws IOException {
		final File barcodes = File.createTempFile("cell_barcodes", ".txt");
		barcodes.deleteOnExit();
		Files.write(barcodes.toPath(), Arrays.asList(CELL_A, CELL_B));
		GenerateSyntheticDoublets g = makeBase();
		g.CELL_BC_FILE = barcodes;
		g.NUMBER_MULTICELL = 1;
		g.MULTIPLICITY = 2;
		Assert.assertEquals(g.doWork(), 0);
		Assert.assertEquals(Files.readAllLines(g.SUMMARY.toPath()).size(), 1);
		Assert.assertEquals(countByTag(g.OUTPUT).values().stream().mapToInt(Integer::intValue).sum(), 8);
	}

	@Test
	public void testRandomMultiCellsFromTopBarcodes() {
		GenerateSyntheticDoublets g = makeBase();
		g.NUM_BARCODES = 2;
		g.NUMBER_MULTICELL = 1;
		g.MULTIPLICITY = 2;
		Assert.assertEquals(g.doWork(), 0);
		Assert.assertEquals(java.nio.file.Files.readAllLines(g.SUMMARY.toPath()).size(), 1);
	}

	@Test(expectedExceptions = IllegalArgumentException.class)
	public void testTooManyMultiCellsRequested() throws IOException {
		final File barcodes = File.createTempFile("cell_barcodes", ".txt");
		barcodes.deleteOnExit();
		Files.write(barcodes.toPath(), Arrays.asList(CELL_A, CELL_B));
		GenerateSyntheticDoublets g = makeBase();
		g.CELL_BC_FILE = barcodes;
		g.NUMBER_MULTICELL = 2;
		g.MULTIPLICITY = 2;
		g.doWork();
	}
}
