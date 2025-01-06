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
package org.broadinstitute.dropseqrna.readtrimming;

import htsjdk.samtools.*;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SequenceBaseEnum;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.SamHeaderUtil;

import java.util.ArrayList;
import java.util.Arrays;

@CommandLineProgramProperties(summary = "Trim homopolymer run from beginning of reads.  If the entire read" +
        " is a homopolymer run, instead of trimming all bases are set to Q3.",
        oneLineSummary = "Trim the given sequence from the beginning of reads",
        programGroup = DropSeq.class)
public class TrimHomopolymerStartingSequence
extends AbstractTrimmerClp {
    public static final byte FULL_TRIM_QUALITY_SCORE = 3;
    private final Log log = Log.getInstance(TrimHomopolymerStartingSequence.class);
    private int readsCompletelyTrimmed = 0;

    @Argument(doc = "How many mismatches are acceptable in the homopolymer run.")
    public int MISMATCHES = 1;

    @Argument(doc = "How many bases of homopolymer qualifies as a run to be trimmed.")
    public int NUM_BASES = 6;

    @Argument(doc = "The base to trim.")
    public SequenceBaseEnum BASE = SequenceBaseEnum.T;
    private byte BASE_BYTE;

    @Override
    protected String[] customCommandLineValidation() {
        final ArrayList<String> list = new ArrayList<>(1);
        if (NUM_BASES < 1) {
            list.add("NUM_BASES must be greater than 0.");
        }
        if (MISMATCHES < 0) {
            list.add("MISMATCHES must be greater than or equal to 0.");
        }
        return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);
    }

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        BASE_BYTE = StringUtil.charToByte(BASE.base);
        if (!SequenceUtil.isValidBase(BASE_BYTE)) {
            throw new IllegalArgumentException(String.format("BASE (%c) must be a valid base, one of ACGT.", BASE));
        }
        final ProgressLogger progress = new ProgressLogger(log);

        SamReader bamReader = SamReaderFactory.makeDefault().open(this.INPUT);
        SAMFileHeader header = bamReader.getFileHeader();
        SamHeaderUtil.addPgRecord(header, this);
        SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT);
        for (SAMRecord r : bamReader) {
            if (shouldTrim(r)) {
                maybeClip(r);
            }
            writer.addAlignment(r);
            progress.record(r);
            this.numReadsTotal++;
        }

        CloserUtil.close(bamReader);

        writer.close();
        log.info("Number of reads trimmed: " + this.readsTrimmed, " total reads: " + this.numReadsTotal);
        if (this.OUTPUT_SUMMARY!=null) writeSummary();

        return 0;
    }

    private void maybeClip(final SAMRecord r) {
        final byte[] readBases = r.getReadBases();
        final int trimPosition = getTrimPosition(readBases);
        if (trimPosition > 0) {
            if (trimPosition == readBases.length) {
                // if the entire read is a homopolymer run, set all qualities to Q3
                byte[] value = new byte[readBases.length];
                Arrays.fill(value, FULL_TRIM_QUALITY_SCORE);
                r.setBaseQualities(value);
                this.readsCompletelyTrimmed++;
            } else {
                final byte[] newReadBases = new byte[readBases.length - trimPosition];
                System.arraycopy(readBases, trimPosition, newReadBases, 0, newReadBases.length);
                r.setReadBases(newReadBases);
                final byte[] readQuals = r.getBaseQualities();
                final byte[] newReadQuals = new byte[readQuals.length - trimPosition];
                System.arraycopy(readQuals, trimPosition, newReadQuals, 0, newReadQuals.length);
                r.setBaseQualities(newReadQuals);
            }
            this.readsTrimmed++;
        }
        this.numBasesTrimmed.increment(trimPosition);
    }

    private int getTrimPosition(final byte[] readBases) {
        if (readBases.length < NUM_BASES) {
            return 0;
        }
        int trimPosition = readBases.length;
        int numMismatches = 0;
        for (int i = 0; i < readBases.length; i++) {
            if (!SequenceUtil.basesEqual(readBases[i], BASE_BYTE)) {
                if (++numMismatches > MISMATCHES) {
                    trimPosition = i;
                    break;
                }
            }
        }
        if (trimPosition < NUM_BASES) {
            return 0;
        }
        if (numMismatches == trimPosition) {
            // handle case in which all that was seen were mismatches.  This shouldn't happen unless MISMATCHES > 1
            return 0;
        }
        return trimPosition;
    }

    private void writeSummary() {
        final TrimMetric metric = new TrimMetric();
        metric.TRIMMED_READS = this.readsTrimmed;
        metric.TOTAL_READS = this.numReadsTotal;
        metric.COMPLETED_TRIMMED_READS = this.readsCompletelyTrimmed;
        MetricsFile<TrimMetric, Integer> mf = new MetricsFile<>();
        mf.addMetric(metric);
        mf.addHistogram(this.numBasesTrimmed);
        mf.write(this.OUTPUT_SUMMARY);
    }

    public static class TrimMetric extends MetricBase {
        public long TRIMMED_READS = 0;
        public long TOTAL_READS = 0;
        public long COMPLETED_TRIMMED_READS = 0;
    }

    public static void main(String[] args) {
        new TrimHomopolymerStartingSequence().instanceMainWithExit(args);
    }
}
