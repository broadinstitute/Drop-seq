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

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SamFiles;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.AbstractSplitBamClp;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import org.broadinstitute.dropseqrna.utils.PairedSamRecordIterator;
import org.broadinstitute.dropseqrna.utils.SamHeaderUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readpairs.ReadPair;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.nio.PicardHtsPath;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;

@CommandLineProgramProperties(
        summary = "Correct edit-distance 1 errors in cell barcodes in scRNA-seq read pairs, in which a region of " +
                "one read of the pair contains the raw cell barcode.  The corrected cell barcode is assigned to the " +
                "read in a tag.  The reads are not altered beyond the addition of tags.",
        oneLineSummary = "Correct edit-distance 1 errors in cell barcodes in scRNA-seq read pairs.",
        programGroup = DropSeq.class
)
public class CorrectScrnaReadPairs
extends CommandLineProgram {

    protected static final Log log = Log.getInstance(CorrectScrnaReadPairs.class);
    protected ProgressLogger progressLogger = new ProgressLogger(log, 10000000);


    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input paired-end SAM or BAM files to " +
            "correct.  They must all have the same sort order", minElements = 1)
    public List<PicardHtsPath> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc="Output SAM or BAM, tagged with corrected cell barcodes.")
    public File OUTPUT;

    @Argument(optional=true, shortName="D", doc="Delete input BAM(s) after splitting them. Default: do not delete input BAM(s).")
    public boolean DELETE_INPUTS = false;

    @Argument(optional = true, shortName = "DI", doc="Delete BAM indices corresponding to input BAMs.  Default: DELETE_INPUTS setting.")
    public Boolean DELETE_INPUT_INDICES;

    @ArgumentCollection
    public CorrectScrnaReadPairsArgumentCollection ARGUMENT_COLLECTION = new CorrectScrnaReadPairsArgumentCollection();

    @Override
    protected int doWork() {
        INPUT = FileListParsingUtils.expandPicardHtsPathList(INPUT);
        final List<Path> inputPaths = PicardHtsPath.toPaths(INPUT);
        inputPaths.stream().forEach(p -> IOUtil.assertFileIsReadable(p));
        IOUtil.assertFileIsWritable(OUTPUT);
        if (DELETE_INPUT_INDICES == null) {
            DELETE_INPUT_INDICES = DELETE_INPUTS;
        }
        // Check that input BAM files can be deleted
        if (DELETE_INPUTS) {
            for (final Path bamFile : inputPaths) {
                IOUtil.assertFileIsWritable(bamFile);
            }
        }

        if (DELETE_INPUT_INDICES) {
            for (final Path bamFile : inputPaths) {
                final Path index = SamFiles.findIndex(bamFile);
                if (index != null && Files.exists(index)) {
                    IOUtil.assertFileIsWritable(index);
                }
            }
        }
        final SamHeaderAndIterator  headerAndIterator = SamFileMergeUtil.mergeInputPaths(inputPaths, true);
        final SAMFileWriter samFileWriter = new SAMFileWriterFactory().setCreateIndex(CREATE_INDEX).
                makeSAMOrBAMWriter(headerAndIterator.header, true, OUTPUT);

        SamHeaderUtil.addPgRecord(headerAndIterator.header, this);
        final BarcodeCorrector barcodeCorrector = new BarcodeCorrector(ARGUMENT_COLLECTION);
        barcodeCorrector.setVERBOSITY(VERBOSITY);
        final PairedSamRecordIterator iterator = new PairedSamRecordIterator(headerAndIterator.iterator);
        for (ReadPair pair: new IterableAdapter<>(iterator)) {
            progressLogger.record(pair.getFirstRead());
            barcodeCorrector.correctReadPair(pair);
            samFileWriter.addAlignment(pair.getFirstRead());
            samFileWriter.addAlignment(pair.getSecondRead());
        }
        CloserUtil.close(headerAndIterator.iterator);
        samFileWriter.close();
        if (ARGUMENT_COLLECTION.METRICS != null) {
            barcodeCorrector.writeMetrics(ARGUMENT_COLLECTION.METRICS, getMetricsFile());
        }

        try {
            if (DELETE_INPUTS) {
                inputPaths.stream().forEach(p -> {
                    try {
                        Files.delete(p);
                    } catch (IOException e) {
                        throw new RuntimeIOException(e);
                    }
                });
            }
            if (DELETE_INPUT_INDICES) {
                for (final Path inputBam : inputPaths) {
                    final Path index = SamFiles.findIndex(inputBam);
                    if (index != null && index.toFile().exists()) {
                        Files.delete(index);

                    }
                }
            }
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }

        return 0;
    }
    /** Stock main method, for testing. */
    public static void main(final String[] args) {
        System.exit(new CorrectScrnaReadPairs().instanceMain(args));
    }
}
