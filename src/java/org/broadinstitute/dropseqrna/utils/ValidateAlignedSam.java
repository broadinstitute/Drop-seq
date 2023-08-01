/*
 * MIT License
 *
 * Copyright 2021 Broad Institute
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

import java.io.File;
import java.util.List;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;


@CommandLineProgramProperties(summary = "Validate BAM file(s)",
        oneLineSummary = "Validate each BAM file in the input to make sure that\n" +
        "- It exists and is readable\n" +
        "- Sort order matches the specified one (default: coordinate)\n" +
        "- There is a non-zero sequences in the sequence dictionary\n" +
        "- Sequence dictionaries for all input BAM(s) are identical\n" +
        "- There are non-zero # of reads\n" +
        "- Optionally, there are non-zero # of paired reads",
        programGroup = DropSeq.class)

public class ValidateAlignedSam extends CommandLineProgram {

    private static final Log log = Log.getInstance(ValidateAlignedSam.class);

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file(s) to validate. This argument can accept wildcards, or a file with the suffix .bam_list that contains the locations of multiple BAM files", minElements = 1)
    public List<File> INPUT_BAM;

    @Argument(shortName = StandardOptionDefinitions.SORT_ORDER_SHORT_NAME)
    public SAMFileHeader.SortOrder SORT_ORDER = SAMFileHeader.SortOrder.coordinate;

    @Argument(doc = "Flag indicating whether to check each input BAM file(s) for the presence of a paired-end read", shortName = "PR")
    public boolean CHECK_CONTAINS_PAIRED_READS = false;

    @Override
    protected int doWork() {
        this.INPUT_BAM = FileListParsingUtils.expandFileList(INPUT_BAM);

        File templateBamFile = null;
        SAMSequenceDictionary templateDictionary = null;

        StringBuilder errorMsg = new StringBuilder();
        SamReaderFactory samReaderFactory = SamReaderFactory.makeDefault();
        for (File inputBam : INPUT_BAM) {
            try {
                IOUtil.assertFileIsReadable(inputBam);
            } catch (Exception ex) {
                appendErrorToErrorMessage(ex.getMessage(), errorMsg);
                continue;
            }

            final SamReader samReader = samReaderFactory.open(inputBam);
            final SAMFileHeader header = samReader.getFileHeader();
            if (header.getSortOrder() != SORT_ORDER) {
                appendErrorToErrorMessage("The expected BAM sort order is " + SORT_ORDER + ", but sort order in " + inputBam.getAbsolutePath() + " is " + header.getSortOrder(), errorMsg);
            }
            if (header.getSequenceDictionary().size() == 0) {
                appendErrorToErrorMessage("Sequence dictionary in " + inputBam.getAbsolutePath() + " is empty", errorMsg);
            } else if (templateDictionary == null) {
                templateBamFile = inputBam;
                templateDictionary = header.getSequenceDictionary();
            } else if (!header.getSequenceDictionary().isSameDictionary(templateDictionary)) {
                appendErrorToErrorMessage("Sequence dictionary in " + inputBam.getAbsolutePath() + " differs from that in " + templateBamFile.getAbsolutePath(), errorMsg);
            }

            boolean bamContainsReads = false;
            boolean bamContainsPairedReads = false;
            for (SAMRecord record : samReader) {
                bamContainsReads = true;
                if (CHECK_CONTAINS_PAIRED_READS) {
                    if (record.getReadPairedFlag()) {
                        bamContainsPairedReads = true;
                        break;
                    }
                } else {
                    break;
                }
            }

            if (!bamContainsReads) {
                appendErrorToErrorMessage("BAM file " + inputBam.getAbsolutePath() + " does not contain any reads", errorMsg);
            }
            if (CHECK_CONTAINS_PAIRED_READS && !bamContainsPairedReads) {
                appendErrorToErrorMessage("BAM file " + inputBam.getAbsolutePath() + " does not contain any paired-end reads", errorMsg);
            }
        }

        if (errorMsg.length() > 0) {
            System.err.println("**********");
            System.err.println(errorMsg);
            System.err.println("**********");
        }

        return errorMsg.length() == 0 ? 0 : 1;
    }

    private void appendErrorToErrorMessage(String error, StringBuilder errorMsg) {
        errorMsg.append(error).append("\n");
    }
}
