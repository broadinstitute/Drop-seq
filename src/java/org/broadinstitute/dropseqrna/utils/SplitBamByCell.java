/*
 * MIT License
 *
 * Copyright 2022 Broad Institute
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
package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.IterableAdapter;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;

import java.io.File;
import java.util.ArrayList;

@CommandLineProgramProperties(
        summary = "Splits input BAM file(s) into NUM_OUTPUTS output BAM files, " +
                "in such a way that all the reads for each cell barcode are in exactly one output BAM file.\n" +
                "Option NUM_OUTPUTS is either supplied explicitly, or computed by dividing the total input BAM file(s) size by the value of TARGET_BAM_SIZE.\n\n" +
                "A typical invocation would look like:\n" +
                "SplitBamByCell I=input.bam OUTPUT=output.__SPLITNUM__.bam OUTPUT_LIST=output.bam_list REPORT=output.split_report.txt NUM_OUTPUTS=25\n\n" +
                "When submitting this command to UGER:\n" +
                "* memory: -m 16G should be sufficient.\n" +
                "* time:   3 minutes/GB of input BAM should be enough.",
        oneLineSummary = "Splits input BAM file(s) by cell barcode",
        programGroup = DropSeq.class
)
public class SplitBamByCell
extends AbstractSplitBamClp {
    @Argument(doc="The tag to examine in order to partition reads.  Set to null to round-robin read (pairs) to output " +
            "files, instead of grouping by cell barcode.")
    public String SPLIT_TAG="XC";

    @Argument(doc =
            "If set to a a value < 1 the program will fail if fewer than this fraction of reads contain SPLIT_TAG." +
                    " If set to a value >= 1, the program will fail if fewer than this many reads contain SPLIT_TAG." +
                    " If set to a negative value then every read must contain SPLIT_TAG.",
            optional = true)
    public double PASSING_READ_THRESHOLD = 0.1;

    @Override
    protected String[] customCommandLineValidation() {
        final ArrayList<String> list = new ArrayList<>();
        if (this.SPLIT_TAG == null && this.OUTPUT_MANIFEST != null) {
            list.add("Cannot produce OUTPUT_MANIFEST if SPLIT_TAG=null.");
        }
        return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);
    }

    protected void splitBAMs () {
        log.info("Splitting BAM files");
        if (SPLIT_TAG == null) {
            // Select writers in a round-robin way, keeping read pairs together
            splitBAMsEqually();
        } else {
            splitBAMsByTag();
        }
    }

    private void splitBAMsEqually() {
        if (headerAndIterator.header.getSortOrder() != SAMFileHeader.SortOrder.queryname) {
            throw new IllegalArgumentException("The input BAM file(s) should be sorted by queryname");
        }


        int writerIdx = -1;
        String lastReadName = null;
        for (SAMRecord r : new IterableAdapter<>(headerAndIterator.iterator)) {
            progressLogger.record(r);

            if (lastReadName == null || !r.getReadName().equals(lastReadName)) {
                writerIdx = (writerIdx + 1) % NUM_OUTPUTS;
                ensureWriter(writerIdx);
            }
            lastReadName = r.getReadName();
            writeRecord(writerIdx, r);
        }
    }

    private void splitBAMsByTag() {
        long numReads = 0;
        long numCellBarcodes = 0;
        final boolean requireCellBarcodes = PASSING_READ_THRESHOLD < 0;
        for (SAMRecord r: new IterableAdapter<>(headerAndIterator.iterator)) {
            progressLogger.record(r);

            final String cellBarcode = r.getStringAttribute(SPLIT_TAG);
            numReads++;
            if (cellBarcode == null) {
                if (requireCellBarcodes) {
                    throw new IllegalArgumentException(
                            "Read " + r.getReadName() + " does not contain the attribute " + SPLIT_TAG
                    );
                }
            } else {
                final int writerIdx = getWriterIdxForCellBarcode(cellBarcode);
                writeRecord(writerIdx, r);
                numCellBarcodes++;
            }
        }

        checkThreshold(numReads, numCellBarcodes);
    }

    private void checkThreshold(long numReads, long numCellBarcodes) {
        log.info(
                String.format(
                        "Processed %d reads of which %d (%.2g) contained the attribute %s and %d (%.2g) did not",
                        numReads,
                        numCellBarcodes,
                        (double) numCellBarcodes / numReads,
                        SPLIT_TAG,
                        numReads - numCellBarcodes,
                        (double) (numReads - numCellBarcodes) / numReads
                )
        );
        if (0 <= PASSING_READ_THRESHOLD && PASSING_READ_THRESHOLD < 1) {
            double pctPassingReads = (double) numCellBarcodes / numReads;
            if (pctPassingReads < PASSING_READ_THRESHOLD) {
                throw new IllegalArgumentException(
                        String.format(
                                "%.2f of reads contained the attribute %s but the threshold is %s",
                                pctPassingReads,
                                SPLIT_TAG,
                                PASSING_READ_THRESHOLD
                        )
                );
            }
        } else if (numCellBarcodes < PASSING_READ_THRESHOLD) {
            throw new IllegalArgumentException(
                    String.format(
                            "%d reads contained the attribute %s but the threshold is %d",
                            numCellBarcodes,
                            SPLIT_TAG,
                            (long) PASSING_READ_THRESHOLD
                    )
            );
        }
    }

    // For backward compatibility with old metrics files
    @SuppressWarnings("unused")
    static class CellBarcodeSplitBamMetric
            extends org.broadinstitute.dropseqrna.utils.CellBarcodeSplitBamMetric {
        public CellBarcodeSplitBamMetric() {
        }

        public CellBarcodeSplitBamMetric(String CELL_BARCODE, int SPLIT_BAM_INDEX, File BAM) {
            super(CELL_BARCODE, SPLIT_BAM_INDEX, BAM);
        }
    }

    // For backward compatibility with old metrics files
    @SuppressWarnings("unused")
    static class SplitBamSummaryMetric
    extends org.broadinstitute.dropseqrna.utils.SplitBamSummaryMetric {
        public SplitBamSummaryMetric(long MEAN, long VARIANCE) {
            super(MEAN, VARIANCE);
        }

        public SplitBamSummaryMetric() {
        }
    }
}
