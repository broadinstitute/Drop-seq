package org.broadinstitute.dropseqrna.utils.readiterators.parallel;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;

public class ParallelBamReaderTest {
    private static final File SMALL_BAM =
            new File("testdata/org/broadinstitute/dropseq/metrics/NucBYReg4Reg.MOUSE.GCTAAGTAAGAT.Elp2.fixed.bam");

    /** Write the records of SMALL_BAM repeated copies times, sorted to match the header. */
    private static Path writeRepeatedBam(final int copies) throws IOException {
        final Path out = Files.createTempFile("ParallelBamReaderTest", ".bam");
        out.toFile().deleteOnExit();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(SMALL_BAM)) {
            final List<SAMRecord> records = new ArrayList<>();
            reader.forEach(records::add);
            try (final SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(reader.getFileHeader(), false, out)) {
                for (int i = 0; i < copies; i++) {
                    records.forEach(r -> writer.addAlignment(r.deepCopy()));
                }
            }
        }
        return out;
    }

    private static List<String> readAll(final SamReader reader) {
        final List<String> result = new ArrayList<>();
        try (final CloseableIterator<SAMRecord> it = reader.iterator()) {
            while (it.hasNext()) {
                result.add(it.next().getSAMString());
            }
        }
        return result;
    }

    private static List<String> readStock(final Path bam) throws IOException {
        try (final SamReader reader = SamReaderFactory.makeDefault().open(bam)) {
            return readAll(reader);
        }
    }

    private static long poolThreadCount() {
        return Thread.getAllStackTraces().keySet().stream()
                .filter(t -> t.getName().startsWith("BamDecodePool-") && t.isAlive()).count();
    }

    @DataProvider(name = "configs")
    public Object[][] configs() {
        // threads, batchSize, depth
        return new Object[][]{
                {1, 1, 0},      // one thread, one record per batch, no read-ahead
                {2, 7, 2},      // batches that do not divide the file evenly
                {4, 1000, 2},   // one batch holds the whole file
                {8, 3, 5},      // more threads than useful, deep read-ahead
        };
    }

    @Test(dataProvider = "configs")
    public void testSameRecordsAsHtsjdk(final int threads, final int batchSize, final int readAhead) throws IOException {
        final Path bam = writeRepeatedBam(30);
        final List<String> expected = readStock(bam);
        Assert.assertEquals(expected.size(), 41 * 30);
        try (final BamDecodePool pool = new BamDecodePool(threads, SamReaderFactory.makeDefault(), batchSize, readAhead, null);
             final SamReader reader = pool.open(bam);
             final SamReader stock = SamReaderFactory.makeDefault().open(bam)) {
            Assert.assertEquals(reader.getFileHeader(), stock.getFileHeader());
            Assert.assertEquals(readAll(reader), expected);
        }
    }

    @Test
    public void testRecordsAreIndependentAndBoundToCallerHeader() throws IOException {
        final Path bam = writeRepeatedBam(20);
        final List<String> expected = readStock(bam);
        try (final BamDecodePool pool = new BamDecodePool(4, SamReaderFactory.makeDefault(), 10, 2, null);
             final SamReader reader = pool.open(bam)) {
            final List<SAMRecord> retained = new ArrayList<>();
            try (final CloseableIterator<SAMRecord> it = reader.iterator()) {
                while (it.hasNext()) {
                    final SAMRecord rec = it.next();
                    Assert.assertSame(rec.getHeader(), reader.getFileHeader());
                    retained.add(rec);
                    // The caller changes its header while workers are decoding.
                    reader.getFileHeader().addComment("comment " + retained.size());
                    reader.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.unsorted);
                }
            }
            Assert.assertEquals(retained.stream().map(SAMRecord::getSAMString).toList(), expected);
        }
    }

    @Test
    public void testOutOfOrderCompletionStillDeliversInOrder() throws Exception {
        final Path bam = writeRepeatedBam(10);
        final List<String> expected = readStock(bam);
        final CountDownLatch batch1Started = new CountDownLatch(1);
        final AtomicInteger completedBeforeBatch0 = new AtomicInteger();
        // Batch 0 may not decode until a later batch has been decoded.
        try (final BamDecodePool pool = new BamDecodePool(4, SamReaderFactory.makeDefault(), 5, 4, batchIndex -> {
            try {
                if (batchIndex == 0) {
                    Assert.assertTrue(batch1Started.await(30, TimeUnit.SECONDS), "batch 1 never started");
                    Thread.sleep(100);
                } else if (batchIndex == 1) {
                    batch1Started.countDown();
                }
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }
        });
             final SamReader reader = pool.open(bam)) {
            Assert.assertEquals(readAll(reader), expected);
        }
    }

    @Test
    public void testManyFilesMergedWithBoundedThreads() throws IOException {
        final List<Path> bams = new ArrayList<>();
        for (int i = 0; i < 50; i++) {
            bams.add(writeRepeatedBam(1 + i % 3));
        }
        for (final boolean maintainSort : new boolean[]{true, false}) {
            final SamHeaderAndIterator serial = SamFileMergeUtil.mergeInputPaths(bams, maintainSort, SamReaderFactory.makeDefault());
            final List<String> expected = new ArrayList<>();
            try (final CloseableIterator<SAMRecord> it = serial.iterator) {
                it.forEachRemaining(r -> expected.add(r.getSAMString()));
            }
            Assert.assertEquals(poolThreadCount(), 0);
            final SamHeaderAndIterator parallel = SamFileMergeUtil.mergeInputPaths(
                    bams, maintainSort, SamReaderFactory.makeDefault(), 3);
            Assert.assertEquals(parallel.header, serial.header);
            final List<String> actual = new ArrayList<>();
            long maxThreads = 0;
            try (final CloseableIterator<SAMRecord> it = parallel.iterator) {
                while (it.hasNext()) {
                    actual.add(it.next().getSAMString());
                    maxThreads = Math.max(maxThreads, poolThreadCount());
                }
            }
            Assert.assertTrue(maxThreads > 0 && maxThreads <= 3, "threads: " + maxThreads);
            Assert.assertEquals(actual.size(), expected.size());
            if (maintainSort) {
                // merge order of equal records is arbitrary in htsjdk, so compare as sorted lists
                Assert.assertEquals(actual.stream().sorted().toList(), expected.stream().sorted().toList());
            } else {
                Assert.assertEquals(actual, expected);
            }
            Assert.assertEquals(poolThreadCount(), 0, "closing the iterator must stop the pool");
        }
    }

    private static List<String> sortKeys(final CloseableIterator<SAMRecord> it) {
        final List<String> keys = new ArrayList<>();
        final htsjdk.samtools.SAMRecordCoordinateComparator comparator = new htsjdk.samtools.SAMRecordCoordinateComparator();
        SAMRecord previous = null;
        while (it.hasNext()) {
            final SAMRecord rec = it.next();
            if (previous != null) {
                Assert.assertTrue(comparator.compare(previous, rec) <= 0,
                        "out of coordinate order: " + previous.getSAMString() + " then " + rec.getSAMString());
            }
            keys.add(rec.getReferenceIndex() + ":" + rec.getAlignmentStart() + ":" + rec.getReadNegativeStrandFlag());
            previous = rec;
        }
        it.close();
        return keys;
    }

    /** A program that needs coordinate order from a sorted merge must get it from the threaded merge too. */
    @Test
    public void testSortedMergeStaysInCoordinateOrder() throws IOException {
        final List<Path> bams = new ArrayList<>();
        for (int i = 0; i < 50; i++) {
            bams.add(writeRepeatedBam(1 + i % 3));
        }
        final List<String> serialKeys =
                sortKeys(SamFileMergeUtil.mergeInputPaths(bams, true, SamReaderFactory.makeDefault()).iterator);
        for (final int numThreads : new int[]{2, 4}) {
            final SamHeaderAndIterator merged =
                    SamFileMergeUtil.mergeInputPaths(bams, true, SamReaderFactory.makeDefault(), numThreads);
            Assert.assertEquals(merged.header.getSortOrder(), SAMFileHeader.SortOrder.coordinate);
            // sortKeys also asserts that consecutive records are never out of order
            Assert.assertEquals(sortKeys(merged.iterator), serialKeys, "numThreads " + numThreads);
        }
    }

    @Test
    public void testAssertSortedWorksOnParallelIterator() throws IOException {
        final Path bam = writeRepeatedBam(10);
        try (final BamDecodePool pool = new BamDecodePool(3);
             final SamReader reader = pool.open(bam)) {
            final CloseableIterator<SAMRecord> it = reader.iterator().assertSorted(SAMFileHeader.SortOrder.coordinate);
            Assert.assertEquals(sortKeys(it).size(), 410);
        }
    }

    /** The pool-wide read-ahead allowance must never stop a merge from getting the next batch it needs. */
    @Test(timeOut = 120000)
    public void testSortedMergeOfManyFilesWithNoReadAhead() throws IOException {
        final List<Path> bams = new ArrayList<>();
        for (int i = 0; i < 50; i++) {
            bams.add(writeRepeatedBam(1 + i % 3));
        }
        final List<String> expected = new ArrayList<>();
        try (final CloseableIterator<SAMRecord> it =
                     SamFileMergeUtil.mergeInputPaths(bams, true, SamReaderFactory.makeDefault()).iterator) {
            it.forEachRemaining(r -> expected.add(r.getSAMString()));
        }
        for (final int readAhead : new int[]{0, 1, 5}) {
            try (final BamDecodePool pool = new BamDecodePool(2, SamReaderFactory.makeDefault(), 3, readAhead, null)) {
                final SamHeaderAndIterator merged = SamFileMergeUtil.mergeInputPaths(bams, true, pool::open);
                final List<String> actual = new ArrayList<>();
                try (final CloseableIterator<SAMRecord> it = merged.iterator) {
                    it.forEachRemaining(r -> actual.add(r.getSAMString()));
                }
                Assert.assertEquals(actual.stream().sorted().toList(), expected.stream().sorted().toList(),
                        "readAhead " + readAhead);
            }
        }
    }

    /** Every read-ahead permit a reader takes must be returned, however the reader ends. */
    @Test
    public void testReadAheadPermitsAreAllReturned() throws IOException {
        final Path bam = writeRepeatedBam(30);
        final int readAhead = 6;
        try (final BamDecodePool pool = new BamDecodePool(3, SamReaderFactory.makeDefault(), 2, readAhead, null)) {
            // one read to the end, one closed early, one never iterated
            try (final SamReader whole = pool.open(bam)) {
                readAll(whole);
            }
            final SamReader early = pool.open(bam);
            final CloseableIterator<SAMRecord> it = early.iterator();
            it.next();
            it.close();
            pool.open(bam).close();
            int acquired = 0;
            while (pool.tryAcquireReadAhead()) {
                acquired++;
            }
            Assert.assertEquals(acquired, readAhead);
        }
    }

    @Test
    public void testHeaderOnlyStartsNoThreads() throws IOException {
        final Path bam = writeRepeatedBam(1);
        try (final BamDecodePool pool = new BamDecodePool(4)) {
            final SamReader reader = pool.open(bam);
            Assert.assertEquals(reader.getFileHeader(), SamReaderFactory.makeDefault().getFileHeader(bam));
            Assert.assertEquals(poolThreadCount(), 0);
        }
    }

    @Test
    public void testCloseVariants() throws IOException {
        final Path bam = writeRepeatedBam(30);
        try (final BamDecodePool pool = new BamDecodePool(4, SamReaderFactory.makeDefault(), 3, 3, null)) {
            // early close, then double close
            final SamReader early = pool.open(bam);
            final CloseableIterator<SAMRecord> it = early.iterator();
            Assert.assertTrue(it.hasNext());
            it.next();
            it.close();
            it.close();
            early.close();
            Assert.assertFalse(it.hasNext());
            // never closed: the pool closes it
            final SamReader neverClosed = pool.open(bam);
            neverClosed.iterator().next();
            // iterator() only once
            final SamReader twice = pool.open(bam);
            twice.iterator();
            Assert.assertThrows(IllegalStateException.class, twice::iterator);
        }
        Assert.assertEquals(poolThreadCount(), 0);
    }

    @Test
    public void testFileOpenedTwiceIsIndependent() throws IOException {
        final Path bam = writeRepeatedBam(5);
        final List<String> expected = readStock(bam);
        try (final BamDecodePool pool = new BamDecodePool(3, SamReaderFactory.makeDefault(), 4, 2, null);
             final SamReader a = pool.open(bam);
             final SamReader b = pool.open(bam)) {
            final CloseableIterator<SAMRecord> ia = a.iterator();
            final CloseableIterator<SAMRecord> ib = b.iterator();
            final List<String> la = new ArrayList<>();
            final List<String> lb = new ArrayList<>();
            // advance a every step and b every other step, then drain b
            while (ia.hasNext()) {
                la.add(ia.next().getSAMString());
                if (la.size() % 2 == 0 && ib.hasNext()) {
                    lb.add(ib.next().getSAMString());
                }
            }
            while (ib.hasNext()) {
                lb.add(ib.next().getSAMString());
            }
            Assert.assertEquals(la, expected);
            Assert.assertEquals(lb, expected);
        }
    }

    @Test
    public void testTruncatedFileFailsInOrder() throws IOException {
        final Path bam = writeRepeatedBam(30);
        final List<String> expected = readStock(bam);
        final byte[] bytes = Files.readAllBytes(bam);
        final Path truncated = Files.createTempFile("ParallelBamReaderTest.truncated", ".bam");
        truncated.toFile().deleteOnExit();
        Files.write(truncated, java.util.Arrays.copyOf(bytes, bytes.length * 6 / 10));
        try (final BamDecodePool pool = new BamDecodePool(4, SamReaderFactory.makeDefault(), 5, 3, null);
             final SamReader reader = pool.open(truncated)) {
            final List<String> got = new ArrayList<>();
            try {
                final CloseableIterator<SAMRecord> it = reader.iterator();
                while (it.hasNext()) {
                    got.add(it.next().getSAMString());
                }
                Assert.fail("Expected a failure reading a truncated BAM; read " + got.size());
            } catch (RuntimeException e) {
                // expected
            }
            Assert.assertTrue(got.size() < expected.size());
            Assert.assertEquals(got, expected.subList(0, got.size()), "records before the failure must be in order");
        }
    }

    /** Records delivered before any failure, and the failure message (null if none). */
    private static class Outcome {
        final List<String> records = new ArrayList<>();
        String error;
    }

    private static Outcome readUntilFailure(final SamReader reader) {
        final Outcome outcome = new Outcome();
        try (final CloseableIterator<SAMRecord> it = reader.iterator()) {
            while (it.hasNext()) {
                outcome.records.add(it.next().getSAMString());
            }
        } catch (RuntimeException | Error e) {   // Error: htsjdk can fail with OutOfMemoryError on a corrupt record size
            outcome.error = e.getClass().getSimpleName() + ": " + e.getMessage();
        }
        return outcome;
    }

    @DataProvider(name = "stringencies")
    public Object[][] stringencies() {
        return new Object[][]{
                {htsjdk.samtools.ValidationStringency.STRICT},
                {htsjdk.samtools.ValidationStringency.LENIENT},
                {htsjdk.samtools.ValidationStringency.SILENT},
        };
    }

    /** An invalid record (mapped, no CIGAR) in the middle of the file is handled as the stock reader handles it. */
    @Test(dataProvider = "stringencies")
    public void testValidationStringencyHonored(final htsjdk.samtools.ValidationStringency stringency) throws IOException {
        final Path bad = Files.createTempFile("ParallelBamReaderTest.invalid", ".bam");
        bad.toFile().deleteOnExit();
        try (final SamReader in = SamReaderFactory.makeDefault().open(SMALL_BAM)) {
            final List<SAMRecord> records = new ArrayList<>();
            in.forEach(records::add);
            try (final SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(in.getFileHeader(), false, bad)) {
                for (int i = 0; i < records.size(); i++) {
                    final SAMRecord rec = records.get(i).deepCopy();
                    if (i == 20) {
                        rec.setReadUnmappedFlag(false);
                        rec.setCigarString("*");
                    }
                    writer.addAlignment(rec);
                }
            }
        }
        final SamReaderFactory factory = SamReaderFactory.makeDefault().validationStringency(stringency);
        final Outcome expected;
        try (final SamReader stock = factory.open(bad)) {
            expected = readUntilFailure(stock);
        }
        if (stringency == htsjdk.samtools.ValidationStringency.STRICT) {
            Assert.assertNotNull(expected.error, "the test file must trigger a validation error");
            Assert.assertTrue(expected.records.size() > 0 && expected.records.size() < 41);
        } else {
            Assert.assertNull(expected.error);
            Assert.assertEquals(expected.records.size(), 41);
        }
        for (final int batchSize : new int[]{1, 5, 1000}) {
            try (final BamDecodePool pool = new BamDecodePool(3, factory, batchSize, 2, null);
                 final SamReader reader = pool.open(bad)) {
                final Outcome actual = readUntilFailure(reader);
                Assert.assertEquals(actual.error, expected.error, "batchSize " + batchSize);
                // htsjdk reads one record ahead and validates it inside next() of the record before it, so it
                // does not hand out the record just before an invalid one.  The parallel reader does.
                final int extra = actual.records.size() - expected.records.size();
                Assert.assertTrue(extra == 0 || (extra == 1 && expected.error != null), "batchSize " + batchSize);
                Assert.assertEquals(actual.records.subList(0, expected.records.size()), expected.records,
                        "batchSize " + batchSize);
            }
        }
    }

    /**
     * A corrupt record size must surface as an error after the records before it, never as a hang (an Error such as
     * OutOfMemoryError on a pool thread would otherwise leave the consumer waiting for a batch that never completes).
     */
    @Test(timeOut = 60000)
    public void testCorruptRecordSizeFailsInsteadOfHanging() throws IOException {
        final Path valid = Files.createTempFile("ParallelBamReaderTest.valid", ".bam");
        valid.toFile().deleteOnExit();
        final List<String> expected = new ArrayList<>();
        try (final SamReader in = SamReaderFactory.makeDefault().open(SMALL_BAM);
             final SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(in.getFileHeader(), true, valid)) {
            final CloseableIterator<SAMRecord> it = in.iterator();
            for (int i = 0; i < 3; i++) {
                final SAMRecord rec = it.next();
                expected.add(rec.getSAMString());
                writer.addAlignment(rec);
            }
        }
        // Inflate the valid BAM, append a record claiming to be Integer.MAX_VALUE bytes long, and deflate it again.
        final java.io.ByteArrayOutputStream inflated = new java.io.ByteArrayOutputStream();
        try (final htsjdk.samtools.util.BlockCompressedInputStream in =
                     new htsjdk.samtools.util.BlockCompressedInputStream(valid)) {
            in.transferTo(inflated);
        }
        inflated.write(new byte[]{(byte) 0xff, (byte) 0xff, (byte) 0xff, 0x7f});
        inflated.write(new byte[100]);
        final Path corrupt = Files.createTempFile("ParallelBamReaderTest.corrupt", ".bam");
        corrupt.toFile().deleteOnExit();
        try (final htsjdk.samtools.util.BlockCompressedOutputStream out =
                     new htsjdk.samtools.util.BlockCompressedOutputStream(corrupt.toFile())) {
            out.write(inflated.toByteArray());
        }
        for (final int batchSize : new int[]{1, 1000}) {
            try (final BamDecodePool pool = new BamDecodePool(3, SamReaderFactory.makeDefault(), batchSize, 2, null);
                 final SamReader reader = pool.open(corrupt)) {
                Assert.assertTrue(reader instanceof ParallelBamReader, "the corrupt file must still be recognized as a BAM");
                final Outcome actual = readUntilFailure(reader);
                Assert.assertNotNull(actual.error, "batchSize " + batchSize);
                // htsjdk reads one record ahead, so it may fail before handing out the last good record.
                Assert.assertTrue(actual.records.size() <= expected.size(), "batchSize " + batchSize);
                Assert.assertEquals(actual.records, expected.subList(0, actual.records.size()),
                        "batchSize " + batchSize + ", error: " + actual.error);
            }
        }
    }

    @Test
    public void testNonBamFallsBackToFactory() throws IOException {
        final Path bam = writeRepeatedBam(2);
        final Path sam = Files.createTempFile("ParallelBamReaderTest", ".sam");
        sam.toFile().deleteOnExit();
        try (final SamReader in = SamReaderFactory.makeDefault().open(bam);
             final SAMFileWriter w = new SAMFileWriterFactory().makeSAMWriter(in.getFileHeader(), true, sam)) {
            in.forEach(w::addAlignment);
        }
        try (final BamDecodePool pool = new BamDecodePool(4);
             final SamReader reader = pool.open(sam)) {
            Assert.assertFalse(reader instanceof ParallelBamReader);
            Assert.assertEquals(readAll(reader).size(), 82);
        }
        Assert.assertEquals(poolThreadCount(), 0);
    }

    @Test
    public void testOneThreadIsStockBehavior() throws IOException {
        final Path bam = writeRepeatedBam(3);
        final SamHeaderAndIterator result = SamFileMergeUtil.mergeInputPaths(
                List.of(bam), true, SamReaderFactory.makeDefault(), 1);
        final List<String> got = new ArrayList<>();
        try (final CloseableIterator<SAMRecord> it = result.iterator) {
            it.forEachRemaining(r -> got.add(r.getSAMString()));
        }
        Assert.assertEquals(got, readStock(bam));
        Assert.assertEquals(poolThreadCount(), 0);
    }
}
