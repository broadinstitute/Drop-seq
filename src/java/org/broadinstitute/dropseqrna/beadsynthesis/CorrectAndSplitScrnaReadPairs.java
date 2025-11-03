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
import org.broadinstitute.barclay.argparser.ArgumentCollection;
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
        summary = "Correct edit-distance 1 errors in cell barcodes in scRNA-seq read pairs, and split into BAMs such " +
                "that all the reads for a cell barcode are in the same BAM.",
        oneLineSummary = "Correct edit-distance 1 errors in cell barcodes in scRNA-seq read pairs.",
        programGroup = DropSeq.class
)
public class CorrectAndSplitScrnaReadPairs
extends AbstractSplitBamClp {
    @ArgumentCollection
    public CorrectScrnaReadPairsArgumentCollection ARGUMENT_COLLECTION = new CorrectScrnaReadPairsArgumentCollection();



    @Override
    protected void splitBAMs() {
        final BarcodeCorrector barcodeCorrector = new BarcodeCorrector(ARGUMENT_COLLECTION);
        barcodeCorrector.setVERBOSITY(VERBOSITY);
        final PairedSamRecordIterator iterator = new PairedSamRecordIterator(headerAndIterator.iterator);
        for (ReadPair pair: new IterableAdapter<>(iterator)) {
            progressLogger.record(pair.getFirstRead());
            final String cellBarcode = barcodeCorrector.correctReadPair(pair);
            final int writerIdx = getWriterIdxForCellBarcode(cellBarcode);
            writeRecord(writerIdx, pair.getFirstRead());
            writeRecord(writerIdx, pair.getSecondRead());
        }

        if (ARGUMENT_COLLECTION.METRICS != null) {
            barcodeCorrector.writeMetrics(ARGUMENT_COLLECTION.METRICS, getMetricsFile());
        }
    }

    @Override
    protected String[] customCommandLineValidation() {
        IOUtil.assertFileIsReadable(ARGUMENT_COLLECTION.ALLOWED_BARCODE_COUNTS);
        final ArrayList<String> list = new ArrayList<>();
        if (ARGUMENT_COLLECTION.BARCODED_READ < 1 || ARGUMENT_COLLECTION.BARCODED_READ > 2) {
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
