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

import htsjdk.samtools.*;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import org.apache.commons.math3.stat.StatUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.dropseqrna.TranscriptomeException;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author skashin
 *
 */

public abstract class AbstractSplitBamClp extends CommandLineProgram {

    protected static final Log log = Log.getInstance(AbstractSplitBamClp.class);
    protected ProgressLogger progressLogger = new ProgressLogger(log, 10000000);


    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM files to analyze.  They must all have the same sort order", minElements = 1)
    public List<File> INPUT;

    @Argument(shortName="N", doc="Number of output files to create", mutex={"TARGET_BAM_SIZE"})
    public Integer NUM_OUTPUTS;

    @Argument(shortName="S", doc="Approximate size of split BAMs to be created. This can be a human-readable number like 500M or 2G", mutex={"NUM_OUTPUTS"})
    public String TARGET_BAM_SIZE;

    @Argument(optional=true, shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Template for output file names.  If OUTPUT_LIST is specified, and OUTPUT is a relative path," +
            " output file paths will be relative to the directory of the OUTPUT_LIST.")
    public File OUTPUT;

    @Argument(optional=true, doc="For each output file, this string in the OUTPUT template will be replaced with an integer.")
    public String OUTPUT_SLUG="__SPLITNUM__";

    @Argument(optional=true, shortName="L", doc="If specified, this file will be created, with NUM_OUTPUTS lines, containing all the output files created.")
    public File OUTPUT_LIST;

    @Argument(optional=true, doc="If specified, this file will be created, containing split BAM files' read count distribution stats.")
    public File REPORT;

    @Argument(optional = true, doc="Write a file with the split index of every cell.  Cannot be used if SPLIT_TAG=null.")
    public File OUTPUT_MANIFEST;

    @Argument(optional=true, shortName="W", doc="Overwrite existing files. Default: fail if any files to be created already exist.")
    public boolean OVERWRITE_EXISTING = false;

    @Argument(optional=true, shortName="D", doc="Delete input BAM(s) after splitting them. Default: do not delete input BAM(s).")
    public boolean DELETE_INPUTS = false;

    @Argument(optional = true, shortName = "DI", doc="Delete BAM indices corresponding to input BAMs.  Default: DELETE_INPUTS setting.")
    public Boolean DELETE_INPUT_INDICES;

    private SAMFileWriterFactory samWriterFactory = null;
    /**
     * History of what cell barcodes go into what writers.  This is updated when a new cell barcode is encountered.
     */
    private final Map<String, Integer> cellBarcodeWriterIdxMap = new HashMap<>();


    protected SamHeaderAndIterator headerAndIterator;
    // Starts out empty, gets expanded to NUM_OUTPUTS writers
    private final List<SAMFileInfo> writerInfoList = new ArrayList<>();

    static final String BAM_EXTENSION = ".bam";
    static final String BAM_LIST_EXTENSION = ".bam_list";
    static final String BAM_REPORT_EXTENSION =".split_bam_report";

    protected void writeRecord(final int writerIdx, final SAMRecord rec) {
        writerInfoList.get(writerIdx).addRecord(rec);
    }

    protected int addWriter() {
        if (writerInfoList.size() >= NUM_OUTPUTS) {
            throw new IllegalStateException("Should not be requesting another writer.");
        }
        int writerIdx = writerInfoList.size();
        writerInfoList.add(createWriterInfo(headerAndIterator.header, writerIdx));
        return writerIdx;
    }

    protected void ensureWriter(final int writerIdx) {
        if (writerIdx == writerInfoList.size()) {
            addWriter();
        } else if (writerIdx > writerInfoList.size()) {
            throw new IllegalStateException("Requesting writer beyond the last one");
        }
        // else nothing to do.
    }

    /**
     * Find the existing writer for the given cell barcode, or create a new writer if we haven't created all of them
     * yet, or pick the existing writer with the fewest reads.
     * @param cellBarcode The barcode to be written.
     * @return index into writerInfoList.
     */
    protected int getWriterIdxForCellBarcode(String cellBarcode) {
        Integer writerIdx = cellBarcodeWriterIdxMap.get(cellBarcode);
        if (writerIdx == null) {
            if (writerInfoList.size() < NUM_OUTPUTS) {
                writerIdx = addWriter();
            } else {
                // find writer with the fewest reads
                Integer minCount = null;
                for (int idx = 0; idx< writerInfoList.size(); idx++) {
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
            throw new TranscriptomeException("Failed to get a writer");
        }
        return writerIdx;
    }

    private enum FileSizeSuffix  {
        k(1024L),
        m(1024L * 1024L),
        g(1024L * 1024L * 1024L),
        t(1024L * 1024L * 1024L * 1024L);

        private final long size;

        FileSizeSuffix(long size) {
            this.size = size;
        }

        public long getSize() {
            return size;
        }
    }

    // public so it can be called from unit test
    @Override
    public int doWork() {
        INPUT = FileListParsingUtils.expandFileList(INPUT);

        if (OUTPUT != null && !OUTPUT.getPath().contains(OUTPUT_SLUG)) {
            throw new IllegalArgumentException(OUTPUT + " does not contain the replacement token " + OUTPUT_SLUG);
        }

        if (OUTPUT == null) {
            setDefaultOutput();
        }

        if (DELETE_INPUT_INDICES == null) {
            DELETE_INPUT_INDICES = DELETE_INPUTS;
        }

        // Check that input BAM files can be deleted
        if (DELETE_INPUTS) {
            for (File bamFile : INPUT) {
                IOUtil.assertFileIsWritable(bamFile);
            }
        }

        if (DELETE_INPUT_INDICES) {
            for (File bamFile : INPUT) {
                final File index = SamFiles.findIndex(bamFile);
                if (index != null && index.exists()) {
                    IOUtil.assertFileIsWritable(index);
                }
            }
        }

        if (TARGET_BAM_SIZE != null) {
            long targetSize = dehumanizeFileSize(TARGET_BAM_SIZE);
            NUM_OUTPUTS = (int) Math.max(1, Math.round(1.0 * getTotalBamSize(INPUT) / targetSize));
        }

        checkOutputOverwrites();

        samWriterFactory = new SAMFileWriterFactory().setCreateIndex(CREATE_INDEX);

        headerAndIterator = SamFileMergeUtil.mergeInputs(INPUT, true);
        SamHeaderUtil.addPgRecord(headerAndIterator.header, this);
        splitBAMs();
        CloserUtil.close(headerAndIterator.iterator);
        for (SAMFileInfo writerInfo : writerInfoList) {
            writerInfo.getWriter().close();
        }
        if (OUTPUT_MANIFEST != null) {
            writeManifest(OUTPUT_MANIFEST);
        }


        writeOutputList();
        writeReport();

        if (DELETE_INPUTS) {
            deleteInputBamFiles();
        }


        return 0;
    }

    private void setDefaultOutput() {
        if (INPUT.size() > 1) {
            throw new IllegalArgumentException("OUTPUT must be specified when INPUT list contains more than one BAM file");
        }

        File bamFile = INPUT.get(0);
        if (!bamFile.getName().endsWith(BAM_EXTENSION)) {
            throw new IllegalArgumentException("Input BAM file " + bamFile.getAbsolutePath() + " does not have the extension " + BAM_EXTENSION);
        }

        String bamRootName = bamFile.getName().replaceAll("\\" + BAM_EXTENSION + "$", "");
        OUTPUT = new File(bamFile.getParent(), bamRootName + "." + OUTPUT_SLUG + BAM_EXTENSION);
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
            actualFileToOpen = FileListParsingUtils.resolveFilePath(OUTPUT_LIST.getAbsoluteFile().getParentFile(), samFile);
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

    protected SAMFileInfo createWriterInfo(final SAMFileHeader header, int splitIdx) {
        final File samFile = getRelativeSplitBamFile(splitIdx);
        final File actualFileToOpen = getActualSplitBamFile(samFile);
        final SAMFileWriter samFileWriter = samWriterFactory.makeSAMOrBAMWriter(header, true, actualFileToOpen);
        return new SAMFileInfo(samFile, samFileWriter, 0);
    }


    /**
     * Split inputs in whatever manner is appropriate
     */
    protected abstract void splitBAMs();

    protected void writeManifest(final File manifest) {
        MetricsFile<CellBarcodeSplitBamMetric, Integer> metricsFile = getMetricsFile();
        for (final Map.Entry<String, Integer> entry : cellBarcodeWriterIdxMap.entrySet()) {
            final Integer fileIndex = entry.getValue();
            metricsFile.addMetric(new CellBarcodeSplitBamMetric(entry.getKey(), fileIndex,
                    writerInfoList.get(fileIndex).samFile));
        }
        final BufferedWriter writer = IOUtil.openFileForBufferedWriting(manifest);
        metricsFile.write(writer);
        try {
            writer.close();
        } catch (IOException e) {
            throw new RuntimeIOException("Exception writing " + manifest.getAbsolutePath(), e);
        }
    }

    private String getOutputNameRoot() {
        String outputName = OUTPUT.getName().replaceAll("\\" + BAM_EXTENSION + "$", "");

        String outputNameRoot;
        int index = outputName.indexOf(OUTPUT_SLUG);
        if (index < 0) {
            outputNameRoot = outputName;
        } else if (index == 0) {
            outputNameRoot = "split";
        } else {
            outputNameRoot = outputName.substring(0, index);
        }
        outputNameRoot = outputNameRoot.replaceAll("\\.+$", "");

        return outputNameRoot;
    }

    private void writeOutputList() {
        final File listFile = (OUTPUT_LIST == null)
            ? new File(OUTPUT.getParent(), getOutputNameRoot() + BAM_LIST_EXTENSION)
            : OUTPUT_LIST;
        final PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(listFile));

        for (SAMFileInfo writerInfo : writerInfoList) {
            out.println(writerInfo.getSamFile().toString());
        }

        out.close();
    }


    private void writeReport() {
        final MetricsFile<SplitBamSummaryMetric, Integer> metricsFile = new MetricsFile<>();

        Histogram<Integer> readCountHistogram = new Histogram<>("SPLIT_INDEX", "READ_COUNT");
        final double[] readCounts = new double[writerInfoList.size()];
        for (int idx=0; idx<writerInfoList.size(); idx++) {
            readCounts[idx] = writerInfoList.get(idx).getReadCount();
            readCountHistogram.increment(idx, readCounts[idx]);
        }
        metricsFile.addMetric(new SplitBamSummaryMetric(StrictMath.round(StatUtils.mean(readCounts)), StrictMath.round(StatUtils.variance(readCounts))));
        metricsFile.addHistogram(readCountHistogram);
        final File reportFile = (REPORT == null)
                ? new File(OUTPUT.getParent(), getOutputNameRoot() + BAM_REPORT_EXTENSION)
                : REPORT;
        metricsFile.write(reportFile);
    }

    private void maybeDeleteIndex(final Path bam) {
        if (DELETE_INPUT_INDICES) {
            final Path index = SamFiles.findIndex(bam);
            if (index != null && index.toFile().exists()) {
                try {
                    Files.delete(index);
                } catch (IOException e) {
                    log.error("Index file " + index.toAbsolutePath() + " could not be deleted.");
                }
            }
        }
    }

    private void maybeDeleteBamFile(final Path bam) {
        maybeDeleteIndex(bam);
        if (DELETE_INPUTS) {
            try {
                Files.delete(bam);
            } catch (IOException e) {
                if (Files.isSymbolicLink(bam)) {
                    // Delete symlink, but don't fail if there is a problem
                    log.error(e, "Symbolic link " + bam.toAbsolutePath() + " could not be deleted");
                } else {
                    throw new RuntimeException("Error deleting " + bam, e);
                }
            }
        }
    }

    private void deleteInputBamFiles() {
        for (File bamFile : INPUT) {
            Path bamPath = bamFile.toPath();
            try {
                while (Files.isSymbolicLink(bamPath)) {
                    Path symlinkPath = bamPath;
                    bamPath = Files.readSymbolicLink(bamPath);
                    maybeDeleteBamFile(symlinkPath);
                }
/*              If input file is in the current directory and does not have a parent, this results in NPE.
                I can't figure out the context in which the code below makes sense, but leaving it here in
                case I'm missing something.
                if (!bamPath.isAbsolute()) {
                    bamPath = bamFile.toPath().getParent().resolve(bamPath);
                }
*/

                maybeDeleteBamFile(bamPath);
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

    protected static class SAMFileInfo {
        private final File samFile;
        private final SAMFileWriter writer;
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

        public void addRecord(final SAMRecord r) {
            writer.addAlignment(r);
            incrementReadCount();
        }
        private void incrementReadCount() {
            readCount++;
        }
    }
}
