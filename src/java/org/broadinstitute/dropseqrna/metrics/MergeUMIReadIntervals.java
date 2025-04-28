/*
 * MIT License
 *
 * Copyright 2023 Broad Institute
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

import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

@CommandLineProgramProperties(summary = "Merges output of GatherUMIReadIntervals.  It is assumed that there is no " +
        "cell barcode overlap between files.",
        oneLineSummary = "Merges output of GatherUMIReadIntervals.",
        programGroup = DropSeq.class)
public class MergeUMIReadIntervals
        extends CommandLineProgram {
    private static final Log log = Log.getInstance(MergeUMIReadIntervals.class);


    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "UMI read intervals files to be merged. " +
            "May be gzipped.",
            minElements = 1)
    public List<File> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Merged UMI read intervals report.  " +
            "May be gzipped.")
    public File OUTPUT;

    @Argument(shortName = "D", doc="Delete input files after merging.")
    public boolean DELETE_INPUTS = false;

    @Override
    protected int doWork() {
        IOUtil.assertFileIsWritable(OUTPUT);
        // Instead of fully parsing the inputs, just use getCurrentLine to output whatever the input is.
        // This will produce nonsensical results if any of the input columns are not in the same order as
        // GatherUMIReadIntervals.writePerUMIStatsHeader defines.
        final MergingIterator<TabbedTextFileWithHeaderParser.Row> it = new MergingIterator<>(new UMIReadIntervalRowComparator(),
                INPUT.stream().map(input -> new TabbedTextFileWithHeaderParser(input).iterator()).collect(Collectors.toList()));
        BufferedWriter out = IOUtil.openFileForBufferedWriting(OUTPUT);
        GatherUMIReadIntervals.writePerUMIStatsHeader(out);
        ProgressLogger prog = new ProgressLogger(log);
        try {
            for (final TabbedTextFileWithHeaderParser.Row row : new IterableAdapter<>(it)) {
                prog.record(row.getField(GatherUMIReadIntervals.UmiReadInterval.CONTIG_HEADER), row.getIntegerField(GatherUMIReadIntervals.UmiReadInterval.POSITION_MIN_HEADER));
                out.write(row.getCurrentLine());
                out.newLine();
            }
            out.close();
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
        CloserUtil.close(it);
        if (DELETE_INPUTS) {
            final List<String> couldNotBeDeleted = INPUT.stream().filter(f -> !f.delete()).map(File::getAbsolutePath).collect(Collectors.toList());
            if (!couldNotBeDeleted.isEmpty()) {
                throw new RuntimeIOException("Unabled to delete files: " + StringUtil.join(", ", couldNotBeDeleted));
            }
        }
        return 0;
    }

    /**
     * It is assumed that input files are in appropriate order so that all that is needed is to order by CBC
     */
    private static class UMIReadIntervalRowComparator
    implements Comparator<TabbedTextFileWithHeaderParser.Row> {
        @Override
        public int compare(TabbedTextFileWithHeaderParser.Row r1, TabbedTextFileWithHeaderParser.Row r2) {
            return r1.getField(GatherUMIReadIntervals.UmiReadInterval.CELL_BARCODE_HEADER).compareTo(r2.getField(GatherUMIReadIntervals.UmiReadInterval.CELL_BARCODE_HEADER));
        }
    }
}
