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
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.BaseRange;
import org.broadinstitute.dropseqrna.utils.SamHeaderUtil;

import java.util.Collections;
import java.util.List;

@CommandLineProgramProperties(summary = "Trim the given base range(s) from all reads of the requested paired-end mode",
        oneLineSummary = "Trim base ranges from reads",
        programGroup = DropSeq.class)
public class ClipReads extends AbstractTrimmerClp {
    private final Log log = Log.getInstance(ClipReads.class);
    @Argument(doc="Base range to remove, separated by a dash.  E.g 1-4.  Can remove multiple ranges by separating them " +
            "with a colon.  For example 1-4:17-22 removes the first 4 bases, then the 17-22 bases.")
    public String BASE_RANGE;

    @Override
    protected int doWork() {
        final ProgressLogger progress = new ProgressLogger(log);
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        SamReader inputSam = SamReaderFactory.makeDefault().open(INPUT);
        SAMFileHeader h= inputSam.getFileHeader();
        SamHeaderUtil.addPgRecord(h, this);
        SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(h, true, OUTPUT);
        List<BaseRange> baseRanges = BaseRange.parseBaseRange(this.BASE_RANGE);
        Collections.sort(baseRanges);
        validateBaseRanges(baseRanges);
        for (SAMRecord r: inputSam) {
            if (shouldTrim(r)) {
                clip(r, baseRanges);
                this.readsTrimmed++;
            }
            writer.addAlignment(r);
            progress.record(r);
            this.numReadsTotal++;
        }
        CloserUtil.close(inputSam);

        writer.close();
        log.info("Total reads: " + this.numReadsTotal, " Reads clipped: " + this.readsTrimmed);

        return 0;
    }

    private void clip(SAMRecord r, List<BaseRange> baseRanges) {
        if (baseRanges.getLast().getEnd() > r.getReadLength()) {
            throw new IllegalArgumentException("Base range end is greater than read length: " + BASE_RANGE + " for read " + r.getReadName());
        }
        baseRanges = BaseRange.invert(baseRanges, r.getReadLength());
        r.setReadBases(BaseRange.getBytesForBaseRange(baseRanges, r.getReadBases()));
        r.setBaseQualities(BaseRange.getBytesForBaseRange(baseRanges, r.getBaseQualities()));
    }

    private void validateBaseRanges(List<BaseRange> baseRanges) {
        int prevEnd = 0;
        for (BaseRange b: baseRanges) {
            if (b.getStart() < 1 || b.getEnd() < 1) {
                throw new IllegalArgumentException("Base ranges must be 1 or greater: " + BASE_RANGE);
            }
            if (b.getStart() > b.getEnd()) {
                throw new IllegalArgumentException("Base range start must be less than or equal to end: " + BASE_RANGE);
            }
            if (b.getStart() <= prevEnd) {
                throw new IllegalArgumentException("Base ranges must not overlap: " + BASE_RANGE);
            }
            prevEnd = b.getEnd();
        }
    }
}
