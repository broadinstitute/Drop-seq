/*
 * MIT License
 *
 * Copyright 2020 Broad Institute
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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.multisample;

import java.io.File;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.metrics.Header;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;

import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import org.broadinstitute.dropseqrna.utils.ReportFileUtil;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.TabbedInputParser;

/**
 * @author skashin
 */
@CommandLineProgramProperties(summary = "Merges doublet detection reports.",
        oneLineSummary = "Merges doublet detection reports.",
        programGroup = DropSeq.class)
public class MergeDoubletAssignments extends CommandLineProgram {

    private static final Log log = Log.getInstance(MergeDoubletAssignments.class);

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input doublet detection reports to be merged.", minElements = 1)
    public List<File> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The merged doublet detection report.")
    public File OUTPUT;


    @Override
    protected int doWork() {
        IOUtil.assertFileIsWritable(OUTPUT);

        INPUT = FileListParsingUtils.expandFileList(INPUT);

        PrintStream writer = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(OUTPUT));
        for (String comment : ReportFileUtil.readComments(INPUT))
            writer.println(comment);
        for (Header header : getDefaultHeaders())
            writer.println(ReportFileUtil.COMMENT_PREFIX + header);

        String header = null;
        TabbedInputParser parser = null;
        final Set<String> cellBarcodes = new HashSet<>();
        try {
            for (File assignmentFile : INPUT) {
                parser = new TabbedInputParser(false, assignmentFile);

                if (!parser.hasNext()) {
                    throw new RuntimeException("Assignments file " + assignmentFile.getAbsolutePath() + " does not contain a header");
                }
                String[] fields = parser.next();
                if (!fields[0].equals("cell")) {
                    throw new RuntimeException("The first field in the header line of the assignments file " + assignmentFile.getAbsolutePath() + " should be 'cell'");
                }
                if (header == null) {
                    header = parser.getCurrentLine();
                    writer.println(header);
                }

                while (parser.hasNext()) {
                    fields = parser.next();
                    String cellBarcode = fields[0];
                    if (cellBarcodes.contains(cellBarcode)) {
                        throw new RuntimeException("Cell barcode " + cellBarcode + " is present in multiple assignment files");
                    }
                    writer.println(parser.getCurrentLine());
                    cellBarcodes.add(cellBarcode);
                }
            }
        } finally {
            parser.close();
            writer.close();
        }

        return 0;
    }

    public static void main(final String[] args) {
        new MergeDoubletAssignments().instanceMainWithExit(args);
    }
}
