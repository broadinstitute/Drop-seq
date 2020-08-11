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
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
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

    @Argument(doc="Number of output files to create", mutex={"TARGET_BAM_SIZE"})
    public Integer NUM_OUTPUTS;

    @Argument(doc="Approximate size of split BAMs to be created. This can be a human-readable number like 500M or 2G")
    public String TARGET_BAM_SIZE;

    @Argument(doc="Template for output file names.  If OUTPUT_LIST is specified, and OUTPUT is a relative path," +
            " output file paths will be relative to the directory of the OUTPUT_LIST.")
    public File OUTPUT;

    @Argument(optional=true, doc="For each output file, this string in the OUTPUT template will be replaced with an integer.")
    public String OUTPUT_SLUG="__SPLITNUM__";

    @Argument(optional=true, doc="If specified, this file will be created, with NUM_OUTPUTS lines, containing all the output files created.")
    public File OUTPUT_LIST;

    @Argument(optional=true, doc="If specified, this file will be created, containing split BAM files' read count distribution stats.")
    public File REPORT;

    @Argument(optional=true, doc="Overwrite existing files. Default: fail if any files to be created already exist.")
    public boolean OVERWRITE_EXISTING = false;

    @Argument(optional=true, doc="Delete input BAM(s) after splitting them. Default: do not delete input BAM(s).")
    public boolean DELETE_INPUTS = false;

    private SAMFileWriterFactory samWriterFactory = null;

    private enum FileSizeSuffix  {
        k(1024L),
        m(1024L * 1024L),
        g(1024L * 1024L * 1024L),
        t(1024L * 1024L * 1024L * 1024L);

        private long size;

        private FileSizeSuffix(long size) {
            this.size = size;
        }

        public long getSize() {
            return size;
        }
    }

    @Override
    protected int doWork() {
        if (!OUTPUT.getPath().contains(OUTPUT_SLUG)) {
            throw new IllegalArgumentException(OUTPUT + " does not contain the replacement token " + OUTPUT_SLUG);
        }

        INPUT = FileListParsingUtils.expandFileList(INPUT);

        // Check that input BAM files can be deleted
        if (DELETE_INPUTS) {
            for (File bamFile : INPUT) {
                IOUtil.assertFileIsWritable(bamFile);
            }
        }

        if (NUM_OUTPUTS == null && TARGET_BAM_SIZE == null || NUM_OUTPUTS != null && TARGET_BAM_SIZE != null) {
            throw new IllegalArgumentException("Error: a value for either NUM_OUTPUTS or TARGET_BAM_SIZE must be specified");
        }

        if (TARGET_BAM_SIZE != null) {
            long targetSize = dehumanizeFileSize(TARGET_BAM_SIZE);
            NUM_OUTPUTS = (int) Math.max(1, Math.round(1.0 * getTotalBamSize(INPUT) / targetSize));
        }

        checkOutputOverwrites();

        samWriterFactory = new SAMFileWriterFactory().setCreateIndex(CREATE_INDEX);

        final List<SAMFileInfo> writerInfoList = new ArrayList<>();
        splitBAMs(writerInfoList);

        if (OUTPUT_LIST != null) {
            writeOutputList(writerInfoList);
        }
        if (REPORT != null) {
            writeReport(writerInfoList);
        }

        if (DELETE_INPUTS) {
            deleteInputBamFiles();
        }

        return 0;
    }

    private File getRelativeSplitBamFile(int splitIdx) {
        final String outputPath = OUTPUT.toString().replace(OUTPUT_SLUG, String.valueOf(splitIdx));
        return new File(outputPath);
    }

    private File getActualSplitBamFile(File samFile) {
        final File actualFileToOpen;
        if (OUTPUT_LIST == null) {
            actualFileToOpen = samFile;
        } else {
            actualFileToOpen = FileListParsingUtils.resolveFilePath(OUTPUT_LIST.getParentFile(), samFile);
        }

        return actualFileToOpen;
    }

    private void checkOutputOverwrites() {
        for (int splitIdx=0; splitIdx<NUM_OUTPUTS; splitIdx++) {
            final File splitBamFile = getActualSplitBamFile(getRelativeSplitBamFile(splitIdx));
            if (splitBamFile.exists()) {
                // Check that this output BAM is not one of the INPUT BAMs
                for (File inputBam : INPUT) {
                    if (splitBamFile.getAbsolutePath().equals(inputBam.getAbsolutePath())) {
                        throw new IllegalArgumentException("Output BAM file " + splitBamFile.getAbsolutePath() + " is the same as input BAM " + inputBam.getAbsolutePath());
                    }
                }

                if (OVERWRITE_EXISTING) {
                    log.warn("BAM file " + splitBamFile.getAbsolutePath() + " exists and will be overwritten");
                } else {
                    throw new IllegalArgumentException("BAM file " + splitBamFile.getAbsolutePath() + " already exists, but OVERWRITE_EXISTING is set to false.");
                }
            }
        }
    }

    private SAMFileInfo createWriterInfo(final SAMFileHeader header, int splitIdx) {
        final File samFile = getRelativeSplitBamFile(splitIdx);
        final File actualFileToOpen = getActualSplitBamFile(samFile);
        final SAMFileWriter samFileWriter = samWriterFactory.makeSAMOrBAMWriter(header, true, actualFileToOpen);
        return new SAMFileInfo(samFile, samFileWriter, 0);
    }

    private void splitBAMsEqually(final List<SAMFileInfo> writerInfoList, final SamHeaderAndIterator headerAndIterator) {
        if (headerAndIterator.header.getSortOrder() != SAMFileHeader.SortOrder.queryname) {
            throw new IllegalArgumentException("The input BAM file(s) should be sorted by queryname");
        }

        ProgressLogger pl = new ProgressLogger(log);

        Integer writerIdx = -1;
        String lastReadName = null;
        for (SAMRecord r : new IterableAdapter<>(headerAndIterator.iterator)) {
            pl.record(r);

            if (lastReadName == null || !r.getReadName().equals(lastReadName)) {
                writerIdx = (writerIdx + 1) % NUM_OUTPUTS;
                if (writerIdx >= writerInfoList.size()) {
                    writerInfoList.add(createWriterInfo(headerAndIterator.header, writerIdx));
                }
            }
            lastReadName = r.getReadName();

            final SAMFileInfo writerInfo = writerInfoList.get(writerIdx);
            writerInfo.getWriter().addAlignment(r);
            writerInfo.incrementReadCount();
        }
    }

    private void splitBAMsByTag(final List<SAMFileInfo> writerInfoList, final SamHeaderAndIterator headerAndIterator) {
        ProgressLogger pl = new ProgressLogger(log);

        final Map<String, Integer> cellBarcodeWriterIdxMap = new HashMap<>();

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
    }

    private void splitBAMs (final List<SAMFileInfo> writerInfoList) {
        log.info("Splitting BAM files");
        final SamHeaderAndIterator headerAndIterator = SamFileMergeUtil.mergeInputs(INPUT, true);
        SamHeaderUtil.addPgRecord(headerAndIterator.header, this);

        if (SPLIT_TAG == null) {
            // Select writers in a round robin way, keeping read pairs together
            splitBAMsEqually(writerInfoList, headerAndIterator);
        } else {
            splitBAMsByTag(writerInfoList, headerAndIterator);
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

    private void deleteInputBamFiles() {
        for (File bamFile : INPUT) {
            Path bamPath = bamFile.toPath();
            try {
                while (Files.isSymbolicLink(bamPath)) {
                    Path symlinkPath = bamPath;
                    bamPath = Files.readSymbolicLink(bamPath);

                    // Delete symlink, but don't fail if there is a problem
                    try {
                        Files.delete(symlinkPath);
                    } catch (IOException var4) {
                        log.error("Symbolic link " + symlinkPath.toAbsolutePath() + " could not be deleted");
                    }
                }
                if (!bamPath.isAbsolute()) {
                    bamPath = bamFile.toPath().getParent().resolve(bamPath);
                }

                Files.delete(bamPath);
            } catch (IOException ex) {
                throw new RuntimeException("Error deleting " + bamPath, ex);
            }
        }
    }

    private long dehumanizeFileSize(String fileSizeString) {
        long conversionFactor = 1;

        char suffix = fileSizeString.charAt(fileSizeString.length() - 1);
        if (Character.isAlphabetic(suffix)) {
            try {
                conversionFactor = FileSizeSuffix.valueOf(String.valueOf(suffix).toLowerCase()).getSize();
                fileSizeString = fileSizeString.substring(0, fileSizeString.length() -1);
            } catch (IllegalArgumentException ex) {
                throw new IllegalArgumentException("Error: " + suffix + " is not a valid file size suffix");
            }
        }

        double size;
        try {
            size = Double.parseDouble(fileSizeString);
        } catch (NumberFormatException ex) {
            throw new IllegalArgumentException("Error: " + fileSizeString + " is not a valid number");
        }

        return (long) (conversionFactor * size);
    }

    private long getTotalBamSize(List<File> inputList) {
        long size = 0;
        try {
            for (File file : inputList) {
                size += Files.size(file.toPath());
            }
        } catch (IOException ex) {
            throw new RuntimeException("Error computing the total BAM file size", ex);
        }

        return size;
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
