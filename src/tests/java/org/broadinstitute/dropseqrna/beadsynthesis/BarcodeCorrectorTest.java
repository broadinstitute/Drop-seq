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
package org.broadinstitute.dropseqrna.beadsynthesis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.broadinstitute.dropseqrna.utils.readpairs.ReadPair;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

public class BarcodeCorrectorTest {
    private static final char[] BASES = {'A', 'C', 'G', 'T'};
    private static final int BARCODE_LENGTH = 8;

    private static String randomBarcode(final Random random, final int length, final char[] alphabet) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < length; ++i) {
            sb.append(alphabet[random.nextInt(alphabet.length)]);
        }
        return sb.toString();
    }

    private static File writeAllowedBarcodeCounts(final Set<String> barcodes, final Random random) {
        final Histogram<String> histogram = new Histogram<>("BARCODE", "COUNT");
        for (final String barcode : barcodes) {
            histogram.increment(barcode, 1 + random.nextInt(1000));
        }
        final MetricsFile<CountBarcodeSequences.CountBarcodeSequenceMetrics, String> metricsFile = new MetricsFile<>();
        metricsFile.addMetric(new CountBarcodeSequences.CountBarcodeSequenceMetrics());
        metricsFile.addHistogram(histogram);
        final File file = TestUtils.getTempReportFile("BarcodeCorrectorTest.", ".barcode_counts");
        metricsFile.write(file);
        return file;
    }

    private static ReadPair makePair(final SAMFileHeader header, final String name, final String barcode, final Random random) {
        final String bases1 = barcode + "TTTTTTTT";
        final SAMRecord r1 = new SAMRecord(header);
        final SAMRecord r2 = new SAMRecord(header);
        for (final SAMRecord r : new SAMRecord[]{r1, r2}) {
            r.setReadName(name);
            r.setReadPairedFlag(true);
            r.setReadUnmappedFlag(true);
            r.setMateUnmappedFlag(true);
        }
        r1.setFirstOfPairFlag(true);
        r2.setSecondOfPairFlag(true);
        r1.setReadString(bases1);
        final byte[] quals = new byte[bases1.length()];
        for (int i = 0; i < quals.length; ++i) {
            quals[i] = (byte) (2 + random.nextInt(38));
        }
        r1.setBaseQualities(quals);
        r2.setReadString("ACGTACGT");
        r2.setBaseQualities(new byte[]{30, 30, 30, 30, 30, 30, 30, 30});
        return new ReadPair(r1, r2);
    }

    /**
     * The packed-barcode fast path must give the same corrected barcodes and metrics as the String-based logic, for
     * exact matches, ED1 corrections (single, multiple and ambiguous), no match, and barcodes with N or lower-case bases.
     */
    @Test
    public void testPackedPathMatchesStringPath() throws IOException {
        final Random random = new Random(7);
        // Few barcodes of length 8 relative to the 65k possible, plus some deliberate ED1 neighbors, so that all
        // correction outcomes occur.
        final Set<String> allowed = new LinkedHashSet<>();
        while (allowed.size() < 3000) {
            allowed.add(randomBarcode(random, BARCODE_LENGTH, BASES));
        }
        final List<String> allowedList = new ArrayList<>(allowed);
        for (int i = 0; i < 300; ++i) {
            final char[] chars = allowedList.get(random.nextInt(allowedList.size())).toCharArray();
            chars[random.nextInt(BARCODE_LENGTH)] = BASES[random.nextInt(4)];
            allowed.add(new String(chars));
        }
        final File counts = writeAllowedBarcodeCounts(allowed, random);

        final CorrectScrnaReadPairsArgumentCollection args = new CorrectScrnaReadPairsArgumentCollection();
        args.BASE_RANGE = "1-" + BARCODE_LENGTH;
        args.ALLOWED_BARCODE_COUNTS = counts;
        args.TAG_BOTH_READS = true;
        args.RAW_BARCODE_TAG = "CR";
        final BarcodeCorrector packed = new BarcodeCorrector(args);
        final BarcodeCorrector reference = new BarcodeCorrector(args);
        reference.disablePackedLookup();

        final SAMFileHeader header = new SAMFileHeader();
        final char[] withN = {'A', 'C', 'G', 'T', 'N', 'a'};
        final java.util.Map<String, Integer> outcomeCounts = new java.util.HashMap<>();
        for (int i = 0; i < 30000; ++i) {
            final String source = allowedList.get(random.nextInt(allowedList.size()));
            char[] chars = source.toCharArray();
            final int kind = random.nextInt(10);
            if (kind < 3) {
                // exact
            } else if (kind < 7) {
                chars[random.nextInt(BARCODE_LENGTH)] = BASES[random.nextInt(4)];
            } else if (kind < 8) {
                chars[random.nextInt(BARCODE_LENGTH)] = withN[random.nextInt(withN.length)];
            } else if (kind < 9) {
                chars[random.nextInt(BARCODE_LENGTH)] = BASES[random.nextInt(4)];
                chars[random.nextInt(BARCODE_LENGTH)] = BASES[random.nextInt(4)];
            } else {
                chars = randomBarcode(random, BARCODE_LENGTH, BASES).toCharArray();
            }
            final String barcode = new String(chars);
            final String name = "read" + i;
            final ReadPair pairA = makePair(header, name, barcode, new Random(i));
            final ReadPair pairB = makePair(header, name, barcode, new Random(i));
            final String fromPacked = packed.correctReadPair(pairA);
            final String fromReference = reference.correctReadPair(pairB);
            Assert.assertEquals(fromPacked, fromReference, barcode);
            Assert.assertEquals(pairA.getFirstRead().getAttribute("XC"), pairB.getFirstRead().getAttribute("XC"));
            Assert.assertEquals(pairA.getSecondRead().getAttribute("CR"), pairB.getSecondRead().getAttribute("CR"));
            outcomeCounts.merge(fromPacked.equals(barcode) ? "unchanged" : "changed", 1, Integer::sum);
        }
        final BarcodeCorrectionMetrics a = packed.getMetrics();
        final BarcodeCorrectionMetrics b = reference.getMetrics();
        Assert.assertEquals(a.NUM_READS_EXACT_MATCH, b.NUM_READS_EXACT_MATCH);
        Assert.assertEquals(a.NUM_READS_CORRECTED_SINGLE_ED1, b.NUM_READS_CORRECTED_SINGLE_ED1);
        Assert.assertEquals(a.NUM_READS_CORRECTED_MULTI_ED1, b.NUM_READS_CORRECTED_MULTI_ED1);
        Assert.assertEquals(a.NUM_READS_UNCORRECTABLE_NO_ED1_BARCODES, b.NUM_READS_UNCORRECTABLE_NO_ED1_BARCODES);
        Assert.assertEquals(a.NUM_READS_UNCORRECTED_AMBIGUOUS, b.NUM_READS_UNCORRECTED_AMBIGUOUS);
        Assert.assertEquals(packed.getNumCandidatesHist(), reference.getNumCandidatesHist());
        // make sure the test exercises every outcome
        Assert.assertTrue(a.NUM_READS_EXACT_MATCH > 0);
        Assert.assertTrue(a.NUM_READS_CORRECTED_SINGLE_ED1 > 0);
        Assert.assertTrue(a.NUM_READS_CORRECTED_MULTI_ED1 > 0);
        Assert.assertTrue(a.NUM_READS_UNCORRECTABLE_NO_ED1_BARCODES > 0);
    }

    /** Allowed barcodes that can't be packed (here an N) turn the fast path off, and correction still works. */
    @Test
    public void testNonAcgtAllowedBarcodeUsesStringPath() {
        final Random random = new Random(3);
        final Set<String> allowed = new LinkedHashSet<>(List.of("ACGTACGN", "ACGTACGT", "TTTTAAAA"));
        final CorrectScrnaReadPairsArgumentCollection args = new CorrectScrnaReadPairsArgumentCollection();
        args.BASE_RANGE = "1-" + BARCODE_LENGTH;
        args.ALLOWED_BARCODE_COUNTS = writeAllowedBarcodeCounts(allowed, random);
        final BarcodeCorrector corrector = new BarcodeCorrector(args);
        final SAMFileHeader header = new SAMFileHeader();
        Assert.assertEquals(corrector.correctReadPair(makePair(header, "a", "ACGTACGN", random)), "ACGTACGN");
        Assert.assertEquals(corrector.correctReadPair(makePair(header, "b", "TTTTAAAC", random)), "TTTTAAAA");
        Assert.assertEquals(corrector.getMetrics().NUM_READS_EXACT_MATCH, 1);
        Assert.assertEquals(corrector.getMetrics().NUM_READS_CORRECTED_SINGLE_ED1, 1);
    }
}
