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

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IterableAdapter;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.AbstractSplitBamClp;
import org.broadinstitute.dropseqrna.utils.PairedSamRecordIterator;
import org.broadinstitute.dropseqrna.utils.readpairs.ReadPair;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.util.ArrayList;

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


    @Override
    protected void splitBAMs() {
        final BarcodeCorrector barcodeCorrector = new BarcodeCorrector(
                ALLOWED_BARCODE_COUNTS,
                BARCODED_READ,
                BASE_RANGE,
                BARCODE_TAG,
                TAG_BOTH_READS,
                LIKELIHOOD_RATIO,
                RAW_BARCODE_TAG,
                BARCODE_QUALS_TAG
        );
        barcodeCorrector.setVERBOSITY(VERBOSITY);
        final PairedSamRecordIterator iterator = new PairedSamRecordIterator(headerAndIterator.iterator);
        for (ReadPair pair: new IterableAdapter<>(iterator)) {
            progressLogger.record(pair.getFirstRead());
            final String cellBarcode = barcodeCorrector.correctReadPair(pair);
            final int writerIdx = getWriterIdxForCellBarcode(cellBarcode);
            writeRecord(writerIdx, pair.getFirstRead());
            writeRecord(writerIdx, pair.getSecondRead());
        }

        if (METRICS != null) {
            final MetricsFile<org.broadinstitute.dropseqrna.beadsynthesis.BarcodeCorrectionMetrics, Integer> metricsFile = getMetricsFile();
            metricsFile.addMetric(barcodeCorrector.getMetrics());
            metricsFile.addHistogram(barcodeCorrector.getNumCandidatesHist());
            metricsFile.write(METRICS);
        }
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

    /** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new CorrectAndSplitScrnaReadPairs().instanceMain(args));
	}

    // To support loading metrics files before this class became top-level
    public static class BarcodeCorrectionMetrics
    extends org.broadinstitute.dropseqrna.beadsynthesis.BarcodeCorrectionMetrics
    {}
}
