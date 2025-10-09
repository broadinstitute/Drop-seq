/*
 * MIT License
 *
 * Copyright 2023 Broad Institute
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
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.AbstractSplitBamClp;
import org.broadinstitute.dropseqrna.utils.BaseRange;
import org.broadinstitute.dropseqrna.utils.PairedSamRecordIterator;
import org.broadinstitute.dropseqrna.utils.StringInterner;
import org.broadinstitute.dropseqrna.utils.readpairs.ReadPair;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.util.*;

@CommandLineProgramProperties(
        summary = "Correct edit-distance 1 errors in cell barcodes in scRNA-seq read pairs.",
        oneLineSummary = "Correct edit-distance 1 errors in cell barcodes in scRNA-seq read pairs.",
        programGroup = DropSeq.class
)
public class CorrectAndSplitScrnaReadPairs
extends AbstractSplitBamClp {
    @Argument(doc = "Which read of each read pair contains the cell barcode [1/2].")
    public int BARCODED_READ = 1;

    @Argument(doc="The region of the barcoded read containing the cell barcode, seperated by a dash.  " +
            "E.g. 1-4.  Can extract multiple ranges by separating them by a colon.  " +
            "For example 1-4:17-22 extracts the first 4 bases, then the 17-22 bases, and glues the sequence together " +
            "into a single cell barcode.")
    public String BASE_RANGE;

    @Argument(doc="Metrics file produced by CountBarcodeSequences that has counts for all the expected cell barcodes " +
            "that are found as exact matches in the input data.")
    public File ALLOWED_BARCODE_COUNTS;

    @Argument(doc="Tag to store the corrected barcode on the non-barcode read.")
    public String BARCODE_TAG = "XC";

    @Argument(doc="If true, assign BARCODE_TAG (also RAW_BARCODE_TAG and BARCODE_QUALS_TAG, if set) to both reads.")
    public boolean TAG_BOTH_READS = false;

    @Argument(shortName = StandardOptionDefinitions.METRICS_FILE_SHORT_NAME, optional = true,
    doc="Various matching and correction metrics")
    public File METRICS;

    @Argument(doc="If more than on allowed barcode matches, (best likelihood)/sum(all likelihoods) " +
            "must be >= this value.")
    public double LIKELIHOOD_RATIO = 0.95;

    @Argument(doc="Store the original barcode sequence in this tag on the non-barcode read.  Default: do not assign this tag.",
    optional = true)
    public String RAW_BARCODE_TAG;
    @Argument(doc="Store the barcode base qualities in this tag on the non-barcode read.  Default: do not assign this tag.",
    optional = true)
    public String BARCODE_QUALS_TAG;

    private Map<String, Double> allowedBarcodeNormalizedOccurences;
    // Don't store many copies of the same allowed barcode in ed1MatchCache
    private final StringInterner allowedBarcodeInterner = new StringInterner();

    private final ResourceLimitedMap<String, List<String>> ed1MatchCache =
            new ResourceLimitedMap<>(
                    1_000_000,
                    new ResourceLimitedMapFunctor<>() {
                        @Override
                        public List<String> makeValue(final String key) {
                            return CorrectAndSplitScrnaReadPairs.this.getEd1Matches(key);
                        }

                        @Override
                        public void finalizeValue(final String key, final List<String> strings) {
                        }
                    }
            );

    private List<BaseRange> baseRanges;

    private final BarcodeCorrectionMetrics metrics = new BarcodeCorrectionMetrics();
    private final Histogram<Integer> numCandidatesHist = new Histogram<>("NUM_ED1_CANDIDATES", "NUM_READS");

    @Override
    protected void splitBAMs() {
        baseRanges = BaseRange.parseBaseRange(this.BASE_RANGE);
        final String baseRangeStr = StringUtil.join(",", baseRanges);
        log.info(String.format("Splitting BAM files based on cell barcode on read %d in range %s",
                BARCODED_READ, baseRangeStr));
        allowedBarcodeNormalizedOccurences = getNormalizedAllowedBarcodes();
        final PairedSamRecordIterator iterator = new PairedSamRecordIterator(headerAndIterator.iterator);
        for (ReadPair pair: new IterableAdapter<>(iterator)) {
            progressLogger.record(pair.getFirstRead());
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
            final int writerIdx = getWriterIdxForCellBarcode(cellBarcode);
            writeRecord(writerIdx, pair.getFirstRead());
            writeRecord(writerIdx, pair.getSecondRead());
        }

        if (METRICS != null) {
            final MetricsFile<BarcodeCorrectionMetrics, Integer> metricsFile = getMetricsFile();
            metricsFile.addMetric(metrics);
            metricsFile.addHistogram(numCandidatesHist);
            metricsFile.write(METRICS);
        }
    }

    private Map<String, Double> getNormalizedAllowedBarcodes() {
        final MetricsFile<CountBarcodeSequences.CountBarcodeSequenceMetrics, String> metricsFile = new MetricsFile<>();
        metricsFile.read(IOUtil.openFileForBufferedReading(ALLOWED_BARCODE_COUNTS));
        Histogram<String> allowedBarcodeHistogram = metricsFile.getHistogram();
        double allowedBarcodeOccurenceCount = allowedBarcodeHistogram.getSumOfValues();
        final Map<String, Double> ret = new HashMap<>(allowedBarcodeHistogram.size());
        for (final String allowedBarcode: allowedBarcodeHistogram.keySet()) {
            ret.put(allowedBarcodeInterner.intern(allowedBarcode), allowedBarcodeHistogram.get(allowedBarcode).getValue()/allowedBarcodeOccurenceCount);
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
                for (final String candidateBarcode: ed1Matches) {
                    double thisLikelihood = computeCandidateBarcodeLikelihood(cellBarcode, candidateBarcode, baseQualities);
                    sumLikelihoods += thisLikelihood;
                    if (thisLikelihood > bestBarcodeLikelihood) {
                        bestBarcodeLikelihood = thisLikelihood;
                        bestBarcode = candidateBarcode;
                    }
                }
                if (bestBarcodeLikelihood/sumLikelihoods >= LIKELIHOOD_RATIO) {
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

    @Override
    protected String[] customCommandLineValidation() {
        IOUtil.assertFileIsReadable(ALLOWED_BARCODE_COUNTS);
        final ArrayList<String> list = new ArrayList<>();
        if (BARCODED_READ < 1 || BARCODED_READ > 2) {
            list.add("BARCODE_READ must be 1 or 2.");
        }
        return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);
    }

    public static class BarcodeCorrectionMetrics extends MetricBase {
        public long NUM_READS_EXACT_MATCH;
        public long NUM_READS_CORRECTED_SINGLE_ED1;
        public long NUM_READS_CORRECTED_MULTI_ED1;
        public long NUM_READS_UNCORRECTABLE_NO_ED1_BARCODES;
        public long NUM_READS_UNCORRECTED_AMBIGUOUS;

        public void merge(final BarcodeCorrectionMetrics that) {
            this.NUM_READS_EXACT_MATCH += that.NUM_READS_EXACT_MATCH;
            this.NUM_READS_CORRECTED_SINGLE_ED1 += that.NUM_READS_CORRECTED_SINGLE_ED1;
            this.NUM_READS_CORRECTED_MULTI_ED1 += that.NUM_READS_CORRECTED_MULTI_ED1;
            this.NUM_READS_UNCORRECTABLE_NO_ED1_BARCODES += that.NUM_READS_UNCORRECTABLE_NO_ED1_BARCODES;
            this.NUM_READS_UNCORRECTED_AMBIGUOUS += that.NUM_READS_UNCORRECTED_AMBIGUOUS;
        }
    }

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new CorrectAndSplitScrnaReadPairs().instanceMain(args));
	}
}
