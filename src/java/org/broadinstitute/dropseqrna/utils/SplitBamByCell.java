/*
 * MIT License
 *
 * Copyright 2019 Broad Institute
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
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.*;
import org.apache.commons.math3.stat.StatUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.TranscriptomeException;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

/**
 *
 * @author skashin
 *
 */

@CommandLineProgramProperties(
        summary = "Splits input BAM file(s) into NUM_OUTPUTS output BAM files, " +
              "in such a way that all the reads for each cell barcode are in exactly one output BAM file",
        oneLineSummary = "Splits input BAM file(s) by cell barcode",
        programGroup = DropSeq.class
)

public class SplitBamByCell extends CommandLineProgram {

    private static final Log log = Log.getInstance(SplitBamByCell.class);

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM files to analyze.  They must all have the same sort order", minElements = 1)
    public List<File> INPUT;

    @Argument(doc="The tag to examine in order to partition reads.")
    public String SPLIT_TAG="XC";

    @Argument(doc="Number of output files to create")
    public Integer NUM_OUTPUTS;

    @Argument(doc="Template for output file names.  If OUTPUT_LIST is specified, and OUTPUT is a relative path," +
            " output file paths will be relative to the directory of the OUTPUT_LIST.")
    public File OUTPUT;

    @Argument(optional=true, doc="For each output file, this string in the OUTPUT template will be replaced with an integer.")
    public String OUTPUT_SLUG="__SPLITNUM__";

    @Argument(optional=true, doc="If specified, this file will be created, with NUM_OUTPUTS lines, containing all the output files created.")
    public File OUTPUT_LIST;

    @Argument(optional=true, doc="If specified, this file will be created, containing split BAM files' read count distribution stats.")
    public File REPORT;

    private SAMFileWriterFactory samWriterFactory = null;

    @Override
    protected int doWork() {
        if (!OUTPUT.getPath().contains(OUTPUT_SLUG)) {
            throw new IllegalArgumentException(OUTPUT + " does not contain the replacement token " + OUTPUT_SLUG);
        }

        samWriterFactory = new SAMFileWriterFactory().setCreateIndex(CREATE_INDEX);

        Map<String, Integer> cellBarcodeWriterIdxMap = new HashMap<>();
        List<SAMFileInfo> writerInfoList = new ArrayList<>();

        splitBAMs(cellBarcodeWriterIdxMap, writerInfoList);

        if (OUTPUT_LIST != null) {
            writeOutputList(writerInfoList);
        }
        if (REPORT != null) {
            writeReport(writerInfoList);
        }

        return 0;
    }

    private SAMFileInfo createWriterInfo(final SAMFileHeader header, int writerIdx) {
        final String outputPath = OUTPUT.toString().replace(OUTPUT_SLUG, String.valueOf(writerIdx));
        final File samFile = new File(outputPath);
        final File actualFileToOpen;
        if (OUTPUT_LIST == null) {
            actualFileToOpen = samFile;
        } else {
            actualFileToOpen = FileListParsingUtils.resolveFilePath(OUTPUT_LIST.getParentFile(), samFile);
        }
        final SAMFileWriter samFileWriter = samWriterFactory.makeSAMOrBAMWriter(header, true, actualFileToOpen);
        return new SAMFileInfo(samFile, samFileWriter, 0);
    }

    private void splitBAMs (final Map<String, Integer> cellBarcodeWriterIdxMap, final List<SAMFileInfo> writerInfoList) {
        log.info("Splitting BAM files");
        final SamHeaderAndIterator headerAndIterator = SamFileMergeUtil.mergeInputs(INPUT, true);
        SamHeaderUtil.addPgRecord(headerAndIterator.header, this);

        ProgressLogger pl = new ProgressLogger(log);
        for (SAMRecord r: new IterableAdapter<>(headerAndIterator.iterator)) {
            pl.record(r);
            final String cellBarcode = r.getStringAttribute(SPLIT_TAG);
            if (cellBarcode == null) {
                throw new IllegalArgumentException("Read " + r.getReadName() + " does not contain the attribute " + SPLIT_TAG);
            }
            Integer writerIdx = cellBarcodeWriterIdxMap.get(cellBarcode);
            if (writerIdx == null) {
                if (writerInfoList.size() < NUM_OUTPUTS) {
                    writerIdx = writerInfoList.size();
                    writerInfoList.add(createWriterInfo(headerAndIterator.header, writerIdx));
                } else {
                    Integer minCount = null;
                    for (int idx=0; idx<writerInfoList.size(); idx++) {
                        int readCount = writerInfoList.get(idx).getReadCount();
                        if (minCount == null || readCount < minCount) {
                            writerIdx = idx;
                            minCount = readCount;
                        }
                    }
                }
                cellBarcodeWriterIdxMap.put(cellBarcode, writerIdx);
            }
            if (writerIdx == null) {
                throw new TranscriptomeException("Failed to get a writer for read " + r.getReadName());
            }
            final SAMFileInfo writerInfo = writerInfoList.get(writerIdx);
            writerInfo.getWriter().addAlignment(r);
            writerInfo.incrementReadCount();
        }
        CloserUtil.close(headerAndIterator.iterator);
        for (SAMFileInfo writerInfo : writerInfoList) {
            writerInfo.getWriter().close();
        }
    }

    private void writeOutputList(final List<SAMFileInfo> writerInfoList) {
        final PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(OUTPUT_LIST));

        for (SAMFileInfo writerInfo : writerInfoList) {
            out.println(writerInfo.getSamFile().toString());
        }

        out.close();
    }

    private void writeReport(final List<SAMFileInfo> writerInfoList) {
        final PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(REPORT));
        out.println("BAM_INDEX" + "\t" + "NUM_READS");

        final double[] readCounts = new double[writerInfoList.size()];
        for (int idx=0; idx<writerInfoList.size(); idx++) {
            readCounts[idx] = writerInfoList.get(idx).getReadCount();
            out.println(idx + "\t" + (int)readCounts[idx]);
        }
        out.println("mean = " + StrictMath.round(StatUtils.mean(readCounts)));
        out.println("variance = " + StrictMath.round(StatUtils.variance(readCounts)));

        out.close();
    }

    private static class SAMFileInfo {
        private File samFile;
        private SAMFileWriter writer;
        int readCount;

        private SAMFileInfo(File samFile, SAMFileWriter writer, int readCount) {
            this.samFile = samFile;
            this.writer = writer;
            this.readCount = readCount;
        }

        public File getSamFile() {
            return samFile;
        }

        public SAMFileWriter getWriter() {
            return writer;
        }

        public int getReadCount() {
            return readCount;
        }

        public void incrementReadCount() {
            readCount++;
        }
    }
}
