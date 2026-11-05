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
package org.broadinstitute.dropseqrna.utils.readiterators.parallel;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Proxy;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.atomic.AtomicInteger;

public class InterleavingSamRecordIteratorTest {
    private static final File SMALL_BAM =
            new File("testdata/org/broadinstitute/dropseq/metrics/NucBYReg4Reg.MOUSE.GCTAAGTAAGAT.Elp2.fixed.bam");
    private static final int SMALL_BAM_RECORDS = 41;

    /** Write the records of SMALL_BAM repeated copies times.  Zero copies gives a header-only file. */
    private static Path writeRepeatedBam(final int copies) throws IOException {
        final Path out = Files.createTempFile("InterleavingSamRecordIteratorTest", ".bam");
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

    private static List<String> readStock(final Path bam) throws IOException {
        final List<String> result = new ArrayList<>();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(bam);
             final CloseableIterator<SAMRecord> it = reader.iterator()) {
            it.forEachRemaining(r -> result.add(r.getSAMString()));
        }
        return result;
    }

    /** Independent round robin: up to chunk records from each list in turn, dropping lists as they run out. */
    private static List<String> expectedRoundRobin(final List<List<String>> inputs, final int chunk) {
        final List<String> result = new ArrayList<>();
        final int[] pos = new int[inputs.size()];
        boolean any = true;
        while (any) {
            any = false;
            for (int i = 0; i < inputs.size(); i++) {
                final List<String> in = inputs.get(i);
                final int end = Math.min(in.size(), pos[i] + chunk);
                if (pos[i] < end) {
                    result.addAll(in.subList(pos[i], end));
                    pos[i] = end;
                    any = true;
                }
            }
        }
        return result;
    }

    private static List<String> drain(final CloseableIterator<SAMRecord> it) {
        final List<String> result = new ArrayList<>();
        while (it.hasNext()) {
            result.add(it.next().getSAMString());
        }
        return result;
    }

    /** A reader that delegates to the real one and counts close() calls. */
    private static SamReader countingReader(final SamReader delegate, final AtomicInteger closes) {
        return (SamReader) Proxy.newProxyInstance(SamReader.class.getClassLoader(), new Class<?>[]{SamReader.class},
                (proxy, method, args) -> {
                    if (method.getName().equals("close")) {
                        closes.incrementAndGet();
                    }
                    try {
                        return method.invoke(delegate, args);
                    } catch (InvocationTargetException e) {
                        throw e.getCause();
                    }
                });
    }

    @Test
    public void testUnequalInputsRoundRobin() throws IOException {
        // 41 * {1, 3, 0, 2} records, chunk 7: the empty input never gets a turn and the short ones drop out early.
        final int[] copies = {1, 3, 0, 2};
        final List<Path> bams = new ArrayList<>();
        final List<List<String>> perFile = new ArrayList<>();
        for (final int c : copies) {
            final Path bam = writeRepeatedBam(c);
            bams.add(bam);
            perFile.add(readStock(bam));
        }
        final List<SamReader> readers = new ArrayList<>();
        for (final Path bam : bams) {
            readers.add(SamReaderFactory.makeDefault().open(bam));
        }
        final List<String> expected = expectedRoundRobin(perFile, 7);
        Assert.assertEquals(expected.size(), SMALL_BAM_RECORDS * Arrays.stream(copies).sum());
        try (final InterleavingSamRecordIterator it = new InterleavingSamRecordIterator(null, readers, 7)) {
            Assert.assertEquals(drain(it), expected);
            Assert.assertFalse(it.hasNext());
        }
    }

    @Test
    public void testHeaderSetOnlyWhenGiven() throws IOException {
        final Path bam = writeRepeatedBam(1);
        final SAMFileHeader outputHeader = new SAMFileHeader();
        try (final SamReader reader = SamReaderFactory.makeDefault().open(bam);
             final InterleavingSamRecordIterator it =
                     new InterleavingSamRecordIterator(outputHeader, List.of(reader), 10)) {
            int n = 0;
            while (it.hasNext()) {
                Assert.assertSame(it.next().getHeader(), outputHeader);
                n++;
            }
            Assert.assertEquals(n, SMALL_BAM_RECORDS);
        }
        try (final SamReader reader = SamReaderFactory.makeDefault().open(bam);
             final InterleavingSamRecordIterator it = new InterleavingSamRecordIterator(null, List.of(reader), 10)) {
            Assert.assertSame(it.next().getHeader(), reader.getFileHeader());
        }
    }

    @Test
    public void testCloseClosesEveryReader() throws IOException {
        final AtomicInteger closes = new AtomicInteger();
        final List<SamReader> readers = new ArrayList<>();
        for (int i = 0; i < 3; i++) {
            readers.add(countingReader(SamReaderFactory.makeDefault().open(writeRepeatedBam(1)), closes));
        }
        final InterleavingSamRecordIterator it = new InterleavingSamRecordIterator(null, readers, 5);
        it.next();
        it.close();
        Assert.assertEquals(closes.get(), 3);
    }

    @Test
    public void testConstructorFailureClosesReaders() throws IOException {
        final AtomicInteger closes = new AtomicInteger();
        final SamReader good = countingReader(SamReaderFactory.makeDefault().open(writeRepeatedBam(1)), closes);
        final SamReader bad = (SamReader) Proxy.newProxyInstance(SamReader.class.getClassLoader(),
                new Class<?>[]{SamReader.class}, (proxy, method, args) -> {
                    if (method.getName().equals("iterator")) {
                        throw new IllegalStateException("cannot start");
                    }
                    if (method.getName().equals("close")) {
                        closes.incrementAndGet();
                    }
                    return null;
                });
        Assert.assertThrows(IllegalStateException.class,
                () -> new InterleavingSamRecordIterator(null, List.of(good, bad), 5));
        Assert.assertEquals(closes.get(), 2);
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testChunkSizeMustBePositive() throws IOException {
        try (final SamReader reader = SamReaderFactory.makeDefault().open(writeRepeatedBam(1))) {
            new InterleavingSamRecordIterator(null, List.of(reader), 0);
        }
    }

    @Test
    public void testFailureSurfacesOnThatInputsTurn() throws IOException {
        final Path good = writeRepeatedBam(30);
        final Path whole = writeRepeatedBam(30);
        final byte[] bytes = Files.readAllBytes(whole);
        final Path truncated = Files.createTempFile("InterleavingSamRecordIteratorTest.truncated", ".bam");
        truncated.toFile().deleteOnExit();
        Files.write(truncated, Arrays.copyOf(bytes, bytes.length * 6 / 10));
        final int batch = 5;
        final List<String> goodRecords = readStock(good);
        try (final BamDecodePool pool = new BamDecodePool(4, SamReaderFactory.makeDefault(), batch, 3, null)) {
            final List<SamReader> readers = List.of(pool.open(good), pool.open(truncated));
            final List<String> got = new ArrayList<>();
            try (final CloseableIterator<SAMRecord> it = pool.interleave(null, readers)) {
                Assert.assertThrows(RuntimeException.class, () -> {
                    while (it.hasNext()) {
                        got.add(it.next().getSAMString());
                    }
                });
            }
            // Whatever came out before the failure is a prefix of the round robin of the two intact inputs.
            final List<String> goodThenTruncatedPrefix = new ArrayList<>();
            final List<String> wholeRecords = readStock(whole);
            final List<String> expectedFull = expectedRoundRobin(List.of(goodRecords, wholeRecords), batch);
            goodThenTruncatedPrefix.addAll(expectedFull.subList(0, got.size()));
            Assert.assertEquals(got, goodThenTruncatedPrefix);
            Assert.assertTrue(got.size() < expectedFull.size());
        }
    }
}
