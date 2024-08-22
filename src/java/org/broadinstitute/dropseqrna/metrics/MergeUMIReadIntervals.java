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
import org.broadinstitute.dropseqrna.metrics.GatherUMIReadIntervals;
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
        final MergingIterator it = new MergingIterator(new UMIReadIntervalComparator(),
                INPUT.stream().map(UMIReadIntervalIterator::new).collect(Collectors.toList()));
        BufferedWriter out = IOUtil.openFileForBufferedWriting(OUTPUT);
        GatherUMIReadIntervals.writePerUMIStatsHeader(out);
        ProgressLogger prog = new ProgressLogger(log);
        for (final GatherUMIReadIntervals.UmiReadInterval interval : new IterableAdapter<GatherUMIReadIntervals.UmiReadInterval>(it)) {
            prog.record(interval.CONTIG, interval.POSITION_MIN);
            GatherUMIReadIntervals.writePerUMIStats(interval, out);
        }
        try {
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
    private static class UMIReadIntervalComparator
    implements Comparator<GatherUMIReadIntervals.UmiReadInterval> {
        @Override
        public int compare(GatherUMIReadIntervals.UmiReadInterval o1, GatherUMIReadIntervals.UmiReadInterval o2) {
            return o1.CELL_BARCODE.compareTo(o2.CELL_BARCODE);
        }
    }

    private static class UMIReadIntervalIterator
    implements CloseableIterator<GatherUMIReadIntervals.UmiReadInterval> {
        private final CloseableIterator<TabbedTextFileWithHeaderParser.Row> iterator;

        public UMIReadIntervalIterator(final File input) {
            iterator = new TabbedTextFileWithHeaderParser(input).iterator();
        }

        @Override
        public void close() {
            iterator.close();
        }

        @Override
        public boolean hasNext() {
            return iterator.hasNext();
        }

        @Override
        public GatherUMIReadIntervals.UmiReadInterval next() {
            return new GatherUMIReadIntervals.UmiReadInterval(iterator.next());
        }
    }
}
