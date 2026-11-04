/*
 * MIT License
 *
 * Copyright 2026 Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
package org.broadinstitute.dropseqrna.vcftools;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.tribble.TribbleException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFEncoder;
import htsjdk.variant.vcf.VCFFileReader;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;

public class ParallelVCFReaderTest {

	// These only bound a hung test; a passing test never waits for them.  They are generous so that a slow or busy
	// machine doesn't turn into a spurious failure.
	private static final long TEST_TIMEOUT_MS = 120_000;
	private static final long HELD_BATCH_WAIT_SECONDS = 60;

	private static final Path DATA_DIR = Path.of("testdata/org/broadinstitute/dropseq/barnyard/digitalallelecounts");
	// 399 records, 23 contigs, 129 samples, bgzipped.
	private static final Path MULTI_CONTIG_VCF = DATA_DIR.resolve("sampleassignment/multisample/TTTGCGCGGAGC:ATTGTTTAGGAG.vcf.gz");
	// 6 records, 2 samples, plain text.
	private static final Path PLAIN_VCF = DATA_DIR.resolve("hek_cells_2snps.vcf");

	/** Every record from the serial htsjdk reader, fully encoded as VCF text, in file order. */
	private static List<String> serialRecords(final Path vcf) {
		try (VCFFileReader reader = new VCFFileReader(vcf, false)) {
			final VCFEncoder encoder = new VCFEncoder(reader.getFileHeader(), true, true);
			final List<String> result = new ArrayList<>();
			for (VariantContext vc : reader) result.add(encoder.encode(vc));
			return result;
		}
	}

	/** Every record from the parallel reader, fully encoded as VCF text, in the order returned. */
	private static List<String> parallelRecords(final ParallelVCFReader reader) {
		final VCFEncoder encoder = new VCFEncoder(reader.getFileHeader(), true, true);
		final List<String> result = new ArrayList<>();
		try (CloseableIterator<VariantContext> it = reader.iterator()) {
			while (it.hasNext()) result.add(encoder.encode(it.next()));
		}
		return result;
	}

	/** Uncompressed BCF, as Picard/htsjdk write it. */
	private Path bcf;

	@BeforeClass
	public void makeBcf() throws IOException {
		bcf = Files.createTempFile("ParallelVCFReaderTest.", ".bcf");
		bcf.toFile().deleteOnExit();
		try (VCFFileReader in = new VCFFileReader(MULTI_CONTIG_VCF, false);
			 VariantContextWriter out = new VariantContextWriterBuilder()
					 .setOutputPath(bcf)
					 .setOutputFileType(VariantContextWriterBuilder.OutputType.BCF)
					 .unsetOption(Options.INDEX_ON_THE_FLY)
					 .build()) {
			out.writeHeader(in.getFileHeader());
			for (VariantContext vc : in) out.add(vc);
		}
	}

	@DataProvider(name = "vcfs")
	public Object[][] vcfs() {
		return new Object[][]{{MULTI_CONTIG_VCF}, {PLAIN_VCF}, {bcf}};
	}

	/**
	 * The same records, with the same content, in the same order as VCFFileReader.
	 * Batch size 10 gives 40 batches (the last one partial) for the 399-record VCF and BCF,
	 * and a single partial batch for the 6-record plain VCF.
	 */
	@Test(dataProvider = "vcfs", timeOut = TEST_TIMEOUT_MS)
	public void testMatchesVCFFileReader(final Path vcf) {
		final int numThreads = 4;
		try (ParallelVCFReader reader = new ParallelVCFReader(vcf, false, numThreads, 10, 2 * numThreads);
			 VCFFileReader serial = new VCFFileReader(vcf, false)) {
			Assert.assertEquals(reader.getFileHeader().getSampleNamesInOrder(), serial.getFileHeader().getSampleNamesInOrder());
			Assert.assertEquals(reader.getFileHeader().getSequenceDictionary(), serial.getFileHeader().getSequenceDictionary());
			final List<String> expected = serialRecords(vcf);
			Assert.assertFalse(expected.isEmpty(), "test fixture has no records");
			Assert.assertEquals(parallelRecords(reader), expected);
		}
	}

	/**
	 * Batches that finish out of order are still returned in file order.
	 * With 2 workers, batch 0 blocks one worker until batch 3 starts.  The other worker must therefore
	 * finish batches 1 and 2 (and start 3) before batch 0 is decoded, so completion is guaranteed to be out of order.
	 */
	@Test(timeOut = TEST_TIMEOUT_MS)
	public void testOutOfOrderCompletionKeepsFileOrder() {
		final CountDownLatch batch3Started = new CountDownLatch(1);
		final AtomicBoolean batch0Waited = new AtomicBoolean(false);
		final int batchSize = 10;
		// Queue depth >= 4 so batches 0-3 can all be in flight while batch 0 is unfinished.
		try (ParallelVCFReader reader = new ParallelVCFReader(MULTI_CONTIG_VCF, false, 2, batchSize, 8, batchIndex -> {
			if (batchIndex == 3) batch3Started.countDown();
			if (batchIndex == 0) {
				try {
					if (!batch3Started.await(HELD_BATCH_WAIT_SECONDS, TimeUnit.SECONDS))
						throw new AssertionError("batch 3 never started while batch 0 was held");
				} catch (InterruptedException e) {
					throw new RuntimeException(e);
				}
				batch0Waited.set(true);
			}
		})) {
			Assert.assertEquals(parallelRecords(reader), serialRecords(MULTI_CONTIG_VCF));
		}
		Assert.assertTrue(batch0Waited.get(), "batch 0 should have been held until batch 3 started");
	}

	/** A VCF with a header and no records. */
	@Test(timeOut = TEST_TIMEOUT_MS)
	public void testHeaderOnlyVcf() throws IOException {
		final Path headerOnly = Files.createTempFile("ParallelVCFReaderTest.headerOnly.", ".vcf");
		headerOnly.toFile().deleteOnExit();
		Files.write(headerOnly, Files.readAllLines(PLAIN_VCF).stream().filter(l -> l.startsWith("#")).toList());
		try (ParallelVCFReader reader = new ParallelVCFReader(headerOnly, false, 4, 10, 8)) {
			Assert.assertEquals(reader.getFileHeader().getSampleNamesInOrder().size(), 2);
			try (CloseableIterator<VariantContext> it = reader.iterator()) {
				Assert.assertFalse(it.hasNext());
			}
		}
	}

	/** A malformed record mid-file throws from the iterator, after the records before it, without hanging. */
	@Test(timeOut = TEST_TIMEOUT_MS)
	public void testMalformedRecordThrows() throws IOException {
		final Path malformed = Files.createTempFile("ParallelVCFReaderTest.malformed.", ".vcf");
		malformed.toFile().deleteOnExit();
		final List<String> lines = new ArrayList<>(Files.readAllLines(PLAIN_VCF));
		int dataLine = 0;
		for (int i = 0; i < lines.size(); i++) {
			if (lines.get(i).startsWith("#")) continue;
			if (dataLine++ == 3) { // corrupt the 4th record's POS
				final String[] fields = lines.get(i).split("\t", -1);
				fields[1] = "notANumber";
				lines.set(i, String.join("\t", fields));
			}
		}
		Files.write(malformed, lines);

		final ParallelVCFReader reader = new ParallelVCFReader(malformed, false, 2, 1, 4);
		final Iterator<VariantContext> it = reader.iterator();
		int returned = 0;
		try {
			while (it.hasNext()) {
				it.next();
				returned++;
			}
			Assert.fail("Expected an exception from the malformed record");
		} catch (RuntimeException expected) {
			Assert.assertEquals(returned, 3, "the three records before the malformed one should be returned");
		} finally {
			reader.close();
		}
		assertNoThreadsAlive(reader);
	}

	/** Closing after reading only part of the file returns promptly and stops every thread. */
	@Test(timeOut = TEST_TIMEOUT_MS)
	public void testCloseEarly() {
		// Small batches and a short queue, so the reader thread is blocked on a full queue when close() is called.
		final ParallelVCFReader reader = new ParallelVCFReader(MULTI_CONTIG_VCF, false, 4, 1, 2);
		final CloseableIterator<VariantContext> it = reader.iterator();
		for (int i = 0; i < 5; i++) it.next();
		it.close();
		Assert.assertFalse(it.hasNext(), "hasNext() after close() should return false, not block");
		assertNoThreadsAlive(reader);
	}

	@DataProvider(name = "oneAndFourThreads")
	public Object[][] oneAndFourThreads() {
		return new Object[][]{{1}, {4}};
	}

	/** requireIndex fails on an unindexed file the same way at any thread count. */
	@Test(dataProvider = "oneAndFourThreads", expectedExceptions = TribbleException.class)
	public void testRequireIndexWithoutIndex(final int numThreads) throws IOException {
		final Path unindexed = Files.createTempFile("ParallelVCFReaderTest.unindexed.", ".vcf");
		unindexed.toFile().deleteOnExit();
		Files.copy(PLAIN_VCF, unindexed, StandardCopyOption.REPLACE_EXISTING);
		new ParallelVCFReader(unindexed, true, numThreads).close();
	}

	/** Like VCFFileReader, the constructor without a requireIndex argument requires an index. */
	@Test(dataProvider = "oneAndFourThreads", expectedExceptions = TribbleException.class)
	public void testDefaultRequiresIndex(final int numThreads) throws IOException {
		final Path unindexed = Files.createTempFile("ParallelVCFReaderTest.unindexedDefault.", ".vcf");
		unindexed.toFile().deleteOnExit();
		Files.copy(PLAIN_VCF, unindexed, StandardCopyOption.REPLACE_EXISTING);
		new ParallelVCFReader(unindexed, numThreads).close();
	}

	/** iterator() may only be called once. */
	@Test(expectedExceptions = IllegalStateException.class)
	public void testIteratorTwiceThrows() {
		try (ParallelVCFReader reader = new ParallelVCFReader(PLAIN_VCF, false, 2)) {
			reader.iterator();
			reader.iterator();
		}
	}

	/** Waits briefly for this reader's threads to exit, then asserts none are still alive. */
	private static void assertNoThreadsAlive(final ParallelVCFReader reader) {
		final String prefix = reader.getThreadNamePrefix();
		final long deadline = System.currentTimeMillis() + 5_000;
		List<String> alive;
		do {
			alive = Thread.getAllStackTraces().keySet().stream()
					.filter(t -> t.getName().startsWith(prefix) && t.isAlive())
					.map(Thread::getName)
					.toList();
			if (alive.isEmpty()) return;
			try {
				Thread.sleep(50);
			} catch (InterruptedException e) {
				throw new RuntimeException(e);
			}
		} while (System.currentTimeMillis() < deadline);
		Assert.fail("Threads still alive after close(): " + alive);
	}
}
