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
package org.broadinstitute.dropseqrna.metrics;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.BaseQualityFilter;
import org.broadinstitute.dropseqrna.utils.TagBamWithReadSequenceExtended;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.File;
import java.util.List;

@CommandLineProgramProperties(summary = "Merges bam_summary.txt files, produced by TagBamWithReadSequenceExtended " +
        "when tagging BAMs with UMIs or CBCs.",
        oneLineSummary = "Merges bam_summary.txt files.",
        programGroup = DropSeq.class)
public class MergeTagBamWithBarcodeSummaries
        extends CommandLineProgram {
    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input rnaseq metrics files to be merged.", minElements = 1)
    public List<File> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The merged rnaseq metrics file.")
    public File OUTPUT;


    @Override
    protected int doWork() {
        IOUtil.assertFileIsWritable(OUTPUT);
        INPUT.forEach(IOUtil::assertFileIsReadable);
        final BaseQualityFilter.FailedBaseMetric outMetric = new BaseQualityFilter.FailedBaseMetric(0);
        for (final File input: INPUT) {
            final TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(input);
            for (final TabbedTextFileWithHeaderParser.Row row: parser) {
                int numFailedBases = row.getIntegerField(TagBamWithReadSequenceExtended.NUM_FAILED_BASES_COLUMN);
                int numBarcodes = row.getIntegerField(TagBamWithReadSequenceExtended.NUM_BARCODES_COLUMN);
                outMetric.addMultipleFailedBase(numFailedBases, numBarcodes);
            }
            CloserUtil.close(parser);
        }
        TagBamWithReadSequenceExtended.writeOutput(outMetric, OUTPUT);
        return 0;
    }
}
