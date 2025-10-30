/*
 * MIT License
 *
 * Copyright 2025 Broad Institute
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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import org.broadinstitute.dropseqrna.utils.BaseRange;
import org.broadinstitute.dropseqrna.utils.StringInterner;
import org.broadinstitute.dropseqrna.utils.readpairs.ReadPair;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Encapsulation of barcode correction part of CorrectAndSplitScrnaReadPairs, so it can be shared with other CLPs.
 */
public class BarcodeCorrector {
    private static final Log log = Log.getInstance(BarcodeCorrector.class);

    // Normalized counts of exact-match barcodes loaded from ALLOWED_BARCODE_COUNTS file
    final private Map<String, Double> allowedBarcodeNormalizedOccurences;
    // Don't store many copies of the same allowed barcode in ed1MatchCache
    private final StringInterner allowedBarcodeInterner = new StringInterner();

    private final ResourceLimitedMap<String, List<String>> ed1MatchCache =
            new ResourceLimitedMap<>(
                    1_000_000,
                    new ResourceLimitedMapFunctor<>() {
                        @Override
                        public List<String> makeValue(final String key) {
                            return getEd1Matches(key);
                        }

                        @Override
                        public void finalizeValue(final String key, final List<String> strings) {
                        }
                    }
            );

    private final List<BaseRange> baseRanges;
    private final int BARCODED_READ;
    private final String BARCODE_TAG;
    private final boolean TAG_BOTH_READS;
    private final double LIKELIHOOD_RATIO;
    private final String RAW_BARCODE_TAG;
    private final String BARCODE_QUALS_TAG;
    private Log.LogLevel VERBOSITY = Log.LogLevel.INFO;
    private final BarcodeCorrectionMetrics metrics = new BarcodeCorrectionMetrics();
    private final Histogram<Integer> numCandidatesHist = new Histogram<>("NUM_ED1_CANDIDATES", "NUM_READS");

    /**
     * See CorrectAndSplitScRnaReadPairs CLP documentation for parameter documentation.
     */
    public BarcodeCorrector(final File ALLOWED_BARCODE_COUNTS,
                            final int BARCODED_READ,
                            final String BASE_RANGE,
                            final String BARCODE_TAG,
                            final boolean TAG_BOTH_READS,
                            final double LIKELIHOOD_RATIO,
                            final String RAW_BARCODE_TAG,
                            final String BARCODE_QUALS_TAG) {
        this.BARCODED_READ = BARCODED_READ;
        this.BARCODE_TAG = BARCODE_TAG;
        this.TAG_BOTH_READS = TAG_BOTH_READS;
        this.LIKELIHOOD_RATIO = LIKELIHOOD_RATIO;
        this.RAW_BARCODE_TAG = RAW_BARCODE_TAG;
        this.BARCODE_QUALS_TAG = BARCODE_QUALS_TAG;
        baseRanges = BaseRange.parseBaseRange(BASE_RANGE);
        String baseRangeStr = StringUtil.join(",", baseRanges);
        log.info(String.format("Correcting cell barcode on read %d in range %s",
                BARCODED_READ, baseRangeStr));
        final MetricsFile<CountBarcodeSequences.CountBarcodeSequenceMetrics, String> metricsFile = new MetricsFile<>();
        metricsFile.read(IOUtil.openFileForBufferedReading(ALLOWED_BARCODE_COUNTS));
        Histogram<String> allowedBarcodeHistogram = metricsFile.getHistogram();
        double allowedBarcodeOccurenceCount = allowedBarcodeHistogram.getSumOfValues();
        allowedBarcodeNormalizedOccurences = new HashMap<>(allowedBarcodeHistogram.size());
        for (final String allowedBarcode: allowedBarcodeHistogram.keySet()) {
            allowedBarcodeNormalizedOccurences.put(allowedBarcodeInterner.intern(allowedBarcode), allowedBarcodeHistogram.get(allowedBarcode).getValue()/allowedBarcodeOccurenceCount);
        }
    }

    /**
     * Tag the read pair according to ctor parameters, and accumulate metrics.
     * @return final CBC, corrected if possible.
     */
    public String correctReadPair(final ReadPair pair) {
        final SAMRecord readWithBarcode = BARCODED_READ == 1? pair.getFirstRead(): pair.getSecondRead();
        final SAMRecord otherRead = BARCODED_READ == 2? pair.getFirstRead(): pair.getSecondRead();
        if (RAW_BARCODE_TAG != null) {
            final String rawBarcode = BaseRange.getSequenceForBaseRange(baseRanges, readWithBarcode.getReadString());
            otherRead.setAttribute(RAW_BARCODE_TAG, rawBarcode);
            if (TAG_BOTH_READS) {
                readWithBarcode.setAttribute(RAW_BARCODE_TAG, rawBarcode);
            }
        }
        if (BARCODE_QUALS_TAG != null) {
            final String barcodeQuals = BaseRange.getSequenceForBaseRange(baseRanges, readWithBarcode.getBaseQualityString());
            otherRead.setAttribute(BARCODE_QUALS_TAG, barcodeQuals);
            if (TAG_BOTH_READS) {
                readWithBarcode.setAttribute(BARCODE_QUALS_TAG, barcodeQuals);
            }
        }
        final String cellBarcode = getCorrectedCellBarcode(readWithBarcode);
        otherRead.setAttribute(BARCODE_TAG, cellBarcode);
        if (TAG_BOTH_READS) {
            readWithBarcode.setAttribute(BARCODE_TAG, cellBarcode);
        }
        return cellBarcode;
    }

    private List<String> getEd1Matches(final String cellBarcode) {
        final List<String> ret = new ArrayList<>();
        final byte[] cellBarcodeBytes = StringUtil.stringToBytes(cellBarcode);
        for (int i = 0; i < cellBarcode.length(); ++i) {
            final byte original = cellBarcodeBytes[i];
            for (byte b: getOtherBases(original)) {
                cellBarcodeBytes[i] = b;
                final String candidate = StringUtil.bytesToString(cellBarcodeBytes);
                if (allowedBarcodeNormalizedOccurences.containsKey(candidate)) {
                    ret.add(allowedBarcodeInterner.intern(candidate));
                }
            }
            cellBarcodeBytes[i] = original;
        }
        return ret;
    }

    private String getCorrectedCellBarcode(final SAMRecord readWithBarcode) {
        final String cellBarcode = BaseRange.getSequenceForBaseRange(baseRanges, readWithBarcode.getReadString());
        if (allowedBarcodeNormalizedOccurences.containsKey(cellBarcode)) {
            // exact match -- no correction needed.
            if (VERBOSITY == Log.LogLevel.DEBUG && metrics.NUM_READS_EXACT_MATCH == 0) {
                log.debug("EXACT_MATCH " + readWithBarcode);
            }
            ++metrics.NUM_READS_EXACT_MATCH;
            return cellBarcode;
        } else {
            List<String> ed1Matches = ed1MatchCache.get(cellBarcode);
            if (ed1Matches.isEmpty()) {
                if (VERBOSITY == Log.LogLevel.DEBUG && metrics.NUM_READS_UNCORRECTABLE_NO_ED1_BARCODES == 0) {
                    log.debug("UNCORRECTABLE_NO_ED1 " + readWithBarcode);
                }
                ++metrics.NUM_READS_UNCORRECTABLE_NO_ED1_BARCODES;
                return cellBarcode;
            } else if (ed1Matches.size() == 1) {
                if (VERBOSITY == Log.LogLevel.DEBUG && metrics.NUM_READS_CORRECTED_SINGLE_ED1 == 0) {
                    log.debug("CORRECTED_SINGLE_ED1 " + readWithBarcode);
                }
                ++metrics.NUM_READS_CORRECTED_SINGLE_ED1;
                numCandidatesHist.increment(1);
                return ed1Matches.getFirst();
            } else {
                String bestBarcode = null;
                double bestBarcodeLikelihood = 0;
                double sumLikelihoods = 0;
                final byte[] baseQualities = BaseRange.getBytesForBaseRange(baseRanges, readWithBarcode.getBaseQualities());
                for (final String candidateBarcode : ed1Matches) {
                    double thisLikelihood = computeCandidateBarcodeLikelihood(cellBarcode, candidateBarcode, baseQualities);
                    sumLikelihoods += thisLikelihood;
                    if (thisLikelihood > bestBarcodeLikelihood) {
                        bestBarcodeLikelihood = thisLikelihood;
                        bestBarcode = candidateBarcode;
                    }
                }
                if (bestBarcodeLikelihood / sumLikelihoods >= LIKELIHOOD_RATIO) {
                    if (VERBOSITY == Log.LogLevel.DEBUG &&
                            (metrics.NUM_READS_CORRECTED_MULTI_ED1 == 0 || ed1Matches.size() >= 4)) {
                        log.debug("CORRECTED_MULTI_ED1 " + readWithBarcode);
                    }
                    ++metrics.NUM_READS_CORRECTED_MULTI_ED1;
                    numCandidatesHist.increment(ed1Matches.size());
                    return bestBarcode;
                } else {
                    if (VERBOSITY == Log.LogLevel.DEBUG && metrics.NUM_READS_UNCORRECTED_AMBIGUOUS == 0) {
                        log.debug("UNCORRECTED_AMBIGUOUS " + readWithBarcode);
                    }
                    ++metrics.NUM_READS_UNCORRECTED_AMBIGUOUS;
                    return cellBarcode;
                }
            }
        }
    }

    private double computeCandidateBarcodeLikelihood(final String uncorrectedBarcode, final String candidateBarcode, final byte[] baseQualities) {
        int i;
        for (i = 0; i < uncorrectedBarcode.length(); ++i) {
            if (uncorrectedBarcode.charAt(i) != candidateBarcode.charAt(i) ) {
                break;
            }
        }
        return QualityUtil.getErrorProbabilityFromPhredScore(baseQualities[i]) * allowedBarcodeNormalizedOccurences.get(candidateBarcode);
    }

    private final byte[] gct = {SequenceUtil.G, SequenceUtil.C, SequenceUtil.T};
    private final byte[] act = {SequenceUtil.A, SequenceUtil.C, SequenceUtil.T};
    private final byte[] agt = {SequenceUtil.A, SequenceUtil.G, SequenceUtil.T};
    private final byte[] acg = {SequenceUtil.A, SequenceUtil.C, SequenceUtil.G};
    private final byte[] acgt = {SequenceUtil.A, SequenceUtil.C, SequenceUtil.G, SequenceUtil.T};
    private byte[] getOtherBases(byte original) {
        return switch (original) {
            case SequenceUtil.A, SequenceUtil.a -> gct;
            case SequenceUtil.C, SequenceUtil.c -> agt;
            case SequenceUtil.G, SequenceUtil.g -> act;
            case SequenceUtil.T, SequenceUtil.t -> acg;
            case SequenceUtil.N, SequenceUtil.n -> acgt;
            default -> throw new IllegalArgumentException(String.format("Unexpected base %d", original));
        };
    }

    public void writeMetrics(final File METRICS,
                             final MetricsFile<BarcodeCorrectionMetrics, Integer> metricsFile) {
        metricsFile.addMetric(this.getMetrics());
        metricsFile.addHistogram(this.getNumCandidatesHist());
        metricsFile.write(METRICS);
    }

    public void setVERBOSITY(Log.LogLevel VERBOSITY) {
        this.VERBOSITY = VERBOSITY;
    }

    /**
     * @return Accumulated metrics of barcodes processed by this object.
     */
    public BarcodeCorrectionMetrics getMetrics() {
        return metrics;
    }

    /**
     * @return histogram with bucket key: number of ED1 matches for a non-exact match barcode;
     * bucket value: number of reads with that number of ED1 matches.
     */
    public Histogram<Integer> getNumCandidatesHist() {
        return numCandidatesHist;
    }
}
