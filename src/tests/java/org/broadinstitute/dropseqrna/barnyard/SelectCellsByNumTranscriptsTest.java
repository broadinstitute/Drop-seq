package org.broadinstitute.dropseqrna.barnyard;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.metrics.MetricsFile;
import org.apache.commons.io.FileUtils;
import org.broadinstitute.dropseqrna.beadsynthesis.GenerateRandomUMIs;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class SelectCellsByNumTranscriptsTest {

	private final File SINGLE_ORGANISM_BAM= new File ("testdata/org/broadinstitute/dropseq/utils/N701_small.bam");
	private final File SINGLE_ORGANISM_BAM_EXPECTED_CB_100_READS=new File ("testdata/org/broadinstitute/dropseq/utils/N701_small.cell_barcodes_100_reads.txt");
	private final File SINGLE_ORGANISM_BAM_EXPECTED_CB_100_TRANSCRIPTS=new File ("testdata/org/broadinstitute/dropseq/utils/N701_small.cell_barcodes_100_transcripts.txt");

	private final File DUAL_ORGANISM_BAM= new File ("testdata/org/broadinstitute/dropseq/utils/human_mouse_smaller.bam");
	private final File DUAL_ORGANISM_BAM_EXPECTED_CB_100_READS=new File ("testdata/org/broadinstitute/dropseq/utils/human_mouse_smaller.cell_barcodes_100_reads.txt");
	private final File DUAL_ORGANISM_BAM_EXPECTED_CB_100_TRANSCRIPTS=new File ("testdata/org/broadinstitute/dropseq/utils/human_mouse_smaller.cell_barcodes_100_transcripts.txt");

	private static SelectCellsByNumTranscripts.Metrics readMetrics(final File file) throws IOException {
		final MetricsFile<SelectCellsByNumTranscripts.Metrics, Integer> metrics = new MetricsFile<>();
		try (Reader in = new FileReader(file)) {
			metrics.read(in);
		}
		return metrics.getMetrics().getFirst();
	}

	private static void assertMetrics(final File file, final int numTranscripts, final int numGenes,
									  final int numTranscriptsSelectedCells, final int numGenesSelectedCells) throws IOException {
		final SelectCellsByNumTranscripts.Metrics m = readMetrics(file);
		Assert.assertEquals(m.NUM_TRANSCRIPTS, numTranscripts);
		Assert.assertEquals(m.NUM_GENES, numGenes);
		Assert.assertEquals(m.NUM_TRANSCRIPTS_SELECTED_CELLS, numTranscriptsSelectedCells);
		Assert.assertEquals(m.NUM_GENES_SELECTED_CELLS, numGenesSelectedCells);
	}

	@Test
	public void testDoWorkSingleOrganism() throws IOException {
		SelectCellsByNumTranscripts s = new SelectCellsByNumTranscripts();
		File outFile = File.createTempFile("SelectCellsByNumTranscripts.", ".cellBarcodes");
		File outFileMinReads = File.createTempFile("SelectCellsByNumTranscripts.", ".minReads");
		outFile.deleteOnExit();
		outFileMinReads.deleteOnExit();

		s.INPUT=Collections.singletonList(this.SINGLE_ORGANISM_BAM);
		s.MIN_TRANSCRIPTS_PER_CELL=100;
		s.OUTPUT=outFile;
		s.MIN_READS_PER_CELL=100;
		s.OUTPUT_INTERIM_CELLS=outFileMinReads;
		s.READ_MQ=0;
		s.ORGANISM=Collections.emptyList();
		s.METRICS = File.createTempFile("SelectCellsByNumTranscripts.", ".metrics");
		s.METRICS.deleteOnExit();

		int success = s.doWork();
		assertMetrics(s.METRICS, 13926, 7259, 12262, 5964);

		try {
			boolean t1 = FileUtils.contentEquals(outFileMinReads, SINGLE_ORGANISM_BAM_EXPECTED_CB_100_READS);
			boolean t2 = FileUtils.contentEquals(outFile, SINGLE_ORGANISM_BAM_EXPECTED_CB_100_TRANSCRIPTS);
			Assert.assertTrue(t1);
			Assert.assertTrue(t2);
		} catch (IOException e) {
			e.printStackTrace();
		}

		Assert.assertEquals(0, success);
	}

	@Test
	public void testDoWorkDualOrganism () throws IOException {

		File outFile = File.createTempFile("SelectCellsByNumTranscripts.", ".cellBarcodes");
		File outFileMinReads = File.createTempFile("SelectCellsByNumTranscripts.", ".minReads");
		outFile.deleteOnExit();
		outFileMinReads.deleteOnExit();

		List<String> organisms = Arrays.asList("HUMAN", "MOUSE");
		SelectCellsByNumTranscripts s = new SelectCellsByNumTranscripts();
		s.INPUT=Collections.singletonList(this.DUAL_ORGANISM_BAM);
		s.ORGANISM=organisms;
		s.MIN_TRANSCRIPTS_PER_CELL=100;
		s.OUTPUT=outFile;
		s.MIN_READS_PER_CELL=null;
		s.OUTPUT_INTERIM_CELLS=outFileMinReads;
		s.READ_MQ=0;
		s.METRICS = File.createTempFile("SelectCellsByNumTranscripts.", ".metrics");
		s.METRICS.deleteOnExit();
		int success = s.doWork();
		assertMetrics(s.METRICS, 94068, 24948, 94011, 24896);

		try {
			boolean t1 = FileUtils.contentEquals(outFileMinReads, DUAL_ORGANISM_BAM_EXPECTED_CB_100_READS);
			boolean t2 = FileUtils.contentEquals(outFile, DUAL_ORGANISM_BAM_EXPECTED_CB_100_TRANSCRIPTS);
			Assert.assertTrue(t1);
			Assert.assertTrue(t2);
		} catch (IOException e) {
			e.printStackTrace();
		}

		Assert.assertEquals(0, success);

	}

	@Test
	public void testReadWriteBarcodes () throws IOException {
		GenerateRandomUMIs g = new GenerateRandomUMIs(0);
		List<String> barcodes = new ArrayList<>();
		for (int i=0; i<20; i++)
			barcodes.add(g.getRandomString(12));
		File outFile = File.createTempFile("SelectCellsByNumTranscripts.", ".cellBarcodes");
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

	@Test
	public void testCustomCommandLineParsing () {
		SelectCellsByNumTranscripts s = new SelectCellsByNumTranscripts();
		List<String> organisms = Arrays.asList("HUMAN", "MOUSE");
		s.ORGANISM=organisms;
		String [] errors = s.customCommandLineValidation();
		Assert.assertNull(errors);

		organisms = Arrays.asList("HUMAN", "HUMAN");
		s.ORGANISM=organisms;
		errors = s.customCommandLineValidation();
		Assert.assertEquals(errors.length, 1);

		organisms = Arrays.asList("HUM"+s.ORGANISM_SEPARATOR+"AN", "MOUSE");
		s.ORGANISM=organisms;
		errors = s.customCommandLineValidation();
		Assert.assertEquals(errors.length, 1);

		s.ORGANISM = Arrays.asList("HUMAN", "MOUSE");
		s.NUM_THREADS = 0;
		errors = s.customCommandLineValidation();
		Assert.assertEquals(errors.length, 1);
	}

	private static final String[] UMIS = {"AAAAAAAA", "CCCCCCCC", "GGGGGGGG"};

	/** Adds one read per UMI for the cell.  The reads are forward strand, so a gene strand of + is sense. */
	private static void addReads(final SAMFileHeader header, final List<SAMRecord> reads, final String cell,
								 final String contig, final String genes, final String strands) {
		for (final String umi : UMIS) {
			final SAMRecord r = new SAMRecord(header);
			r.setReadName(cell + "_" + contig + "_" + umi);
			r.setReferenceName(contig);
			r.setAlignmentStart(100);
			r.setCigarString("10M");
			r.setReadString("ACGTACGTAC");
			r.setBaseQualityString("IIIIIIIIII");
			r.setMappingQuality(255);
			r.setAttribute("XC", cell);
			r.setAttribute("XM", umi);
			r.setAttribute("gn", genes);
			r.setAttribute("gs", strands);
			r.setAttribute("gf", genes.contains(",") ? "CODING,CODING" : "CODING");
			reads.add(r);
		}
	}

	/**
	 * Reads that overlap more than one gene have a comma separated gene tag.  The organism prefix has to be on every
	 * gene in the list, or a gene that is not first is not recognised as belonging to the organism.
	 */
	@Test(dataProvider = "numThreads")
	public void testMultiGeneTagDualOrganism(final int numThreads) throws IOException {
		final SAMFileHeader header = new SAMFileHeader();
		header.addSequence(new SAMSequenceRecord("HUMAN_1", 1000));
		header.addSequence(new SAMSequenceRecord("MOUSE_1", 1000));
		final List<SAMRecord> reads = new ArrayList<>();
		// Only B is on the read's strand, so the read resolves to the second gene in the list.
		addReads(header, reads, "AAAAAAAAAAAA", "HUMAN_1", "A,B", "-,+");
		addReads(header, reads, "CCCCCCCCCCCC", "MOUSE_1", "C,D", "-,+");
		// Both entries are the same gene, so the read resolves to that gene.
		addReads(header, reads, "GGGGGGGGGGGG", "HUMAN_1", "E,E", "+,+");

		final File bam = File.createTempFile("SelectCellsByNumTranscripts.", ".bam");
		bam.deleteOnExit();
		try (SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header, true, bam)) {
			reads.forEach(writer::addAlignment);
		}

		final File outFile = File.createTempFile("SelectCellsByNumTranscripts.", ".cellBarcodes");
		outFile.deleteOnExit();
		final SelectCellsByNumTranscripts s = new SelectCellsByNumTranscripts();
		s.INPUT = Collections.singletonList(bam);
		s.ORGANISM = Arrays.asList("HUMAN", "MOUSE");
		s.MIN_TRANSCRIPTS_PER_CELL = 3;
		s.OUTPUT = outFile;
		s.READ_MQ = 0;
		s.NUM_THREADS = numThreads;
		s.METRICS = File.createTempFile("SelectCellsByNumTranscripts.", ".metrics");
		s.METRICS.deleteOnExit();
		Assert.assertEquals(s.doWork(), 0);

		final List<String> cells = new ArrayList<>(SelectCellsByNumTranscripts.readBarcodes(outFile));
		Collections.sort(cells);
		Assert.assertEquals(cells, Arrays.asList("AAAAAAAAAAAA", "CCCCCCCCCCCC", "GGGGGGGGGGGG"));

		assertMetrics(s.METRICS, 9, 3, 9, 3);
	}

	/** No read has a gene, so no UMI group is produced.  The tool should write an empty cell list, not fail. */
	@Test
	public void testNoGenesDualOrganism() throws IOException {
		final SAMFileHeader header = new SAMFileHeader();
		header.addSequence(new SAMSequenceRecord("HUMAN_1", 1000));
		header.addSequence(new SAMSequenceRecord("MOUSE_1", 1000));
		final List<SAMRecord> reads = new ArrayList<>();
		addReads(header, reads, "AAAAAAAAAAAA", "HUMAN_1", "A", "+");
		reads.forEach(r -> r.setAttribute("gn", null));
		final File bam = File.createTempFile("SelectCellsByNumTranscripts.", ".bam");
		bam.deleteOnExit();
		try (SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header, true, bam)) {
			reads.forEach(writer::addAlignment);
		}
		final File outFile = File.createTempFile("SelectCellsByNumTranscripts.", ".cellBarcodes");
		outFile.deleteOnExit();
		final SelectCellsByNumTranscripts s = new SelectCellsByNumTranscripts();
		s.INPUT = Collections.singletonList(bam);
		s.ORGANISM = Arrays.asList("HUMAN", "MOUSE");
		s.MIN_TRANSCRIPTS_PER_CELL = 3;
		s.OUTPUT = outFile;
		s.READ_MQ = 0;
		Assert.assertEquals(s.doWork(), 0);
		Assert.assertTrue(SelectCellsByNumTranscripts.readBarcodes(outFile).isEmpty());
	}

	/** Cells with the same number of transcripts are in barcode order, so the output does not depend on hash order. */
	@Test
	public void testTiedCellsAreInBarcodeOrder() throws IOException {
		final SAMFileHeader header = new SAMFileHeader();
		header.addSequence(new SAMSequenceRecord("HUMAN_1", 1000));
		header.addSequence(new SAMSequenceRecord("HUMAN_2", 1000));
		final List<SAMRecord> reads = new ArrayList<>();
		final java.util.Random random = new java.util.Random(7);
		final java.util.TreeSet<String> tiedCells = new java.util.TreeSet<>();
		while (tiedCells.size() < 30) {
			final StringBuilder cell = new StringBuilder();
			for (int i = 0; i < 12; i++)
				cell.append("ACGT".charAt(random.nextInt(4)));
			// the cell with more transcripts has the last barcode in order, and must still be first
			if (!cell.toString().equals("TTTTTTTTTTTT"))
				tiedCells.add(cell.toString());
		}
		// three UMIs in one gene is three transcripts
		for (final String cell : tiedCells)
			addReads(header, reads, cell, "HUMAN_1", "A", "+");
		addReads(header, reads, "TTTTTTTTTTTT", "HUMAN_1", "A", "+");
		// a different contig gives the reads different names, so they are not alignments of the same reads
		addReads(header, reads, "TTTTTTTTTTTT", "HUMAN_2", "B", "+");

		final File bam = File.createTempFile("SelectCellsByNumTranscripts.", ".bam");
		bam.deleteOnExit();
		try (SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(header, true, bam)) {
			reads.forEach(writer::addAlignment);
		}
		final File outFile = File.createTempFile("SelectCellsByNumTranscripts.", ".cellBarcodes");
		outFile.deleteOnExit();
		final SelectCellsByNumTranscripts s = new SelectCellsByNumTranscripts();
		s.INPUT = Collections.singletonList(bam);
		s.MIN_TRANSCRIPTS_PER_CELL = 3;
		s.OUTPUT = outFile;
		s.READ_MQ = 0;
		Assert.assertEquals(s.doWork(), 0);

		final List<String> expected = new ArrayList<>();
		expected.add("TTTTTTTTTTTT");
		expected.addAll(tiedCells);
		Assert.assertEquals(SelectCellsByNumTranscripts.readBarcodes(outFile), expected);
	}

	@DataProvider(name = "numThreads")
	public Object[][] numThreads() {
		return new Object[][]{{1}, {2}};
	}
}
