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
package org.broadinstitute.dropseqrna.metrics;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IterableAdapter;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.PositionalArguments;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.TabbedInputParser;

import java.io.File;
import java.io.FileFilter;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.LongAdder;
import java.util.stream.Collectors;

@CommandLineProgramProperties(summary = "Tally frequency of sample indices in barcode.txt.gz files that do not match" +
        " any of the expected sample indices.  This can be useful in diagnosing sample index mixups.",
        oneLineSummary = "Tally frequency of unmatched sample indices.",
        programGroup = DropSeq.class)
public class CountUnmatchedSampleIndices
        extends CommandLineProgram {
    private static final Log LOG = Log.getInstance(CountUnmatchedSampleIndices.class);

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which metrics will be written")
    public File OUTPUT;

    @Argument(doc="This many of the most frequent sample indices will be written to the output file. " +
            "Set to null to remove this limit.",
            optional = true)
    public Integer MAX_OUTPUT = 100;

    @Argument(doc="Only emit sample indices that appear at least this many times.  Default: no threshold.",
            optional = true)
    public Integer MIN_COUNT;

    @Argument(doc="Use this many background threads.", minValue = 1)
    public int NUM_THREADS = 1;

    @PositionalArguments(minElements = 1,
            doc="barcode files (as produced by ExtractIlluminaBarcodes) to be read, or barcode directories.")
    public List<File> BARCODE_FILES;

    @Override
    protected int doWork() {
        final ConcurrentHashMap<String, LongAdder> tally = new ConcurrentHashMap<String, LongAdder>();

        final List<File> allBarcodeFiles = new ArrayList<>(BARCODE_FILES.size());
        final BarcodeFileFilter bff = new BarcodeFileFilter();
        for (final File f : BARCODE_FILES) {
            if (f.isDirectory()) {
                allBarcodeFiles.addAll(Arrays.asList(f.listFiles(bff)));
            } else {
                allBarcodeFiles.add(f);
            }
        }

        if (allBarcodeFiles.isEmpty()) {
            throw new RuntimeException("No barcode files found");
        }
        NUM_THREADS = Math.min(NUM_THREADS, allBarcodeFiles.size());

        final Queue<File> barcodeFiles;

        if (NUM_THREADS > 1) {
            barcodeFiles = new ArrayBlockingQueue<>(allBarcodeFiles.size(), false, allBarcodeFiles);
        } else {
            barcodeFiles = new LinkedList<>(allBarcodeFiles);
        }
        final Worker[] workers = new Worker[NUM_THREADS];
        for (int i = 0; i < NUM_THREADS; ++i) {
            workers[i] = new Worker(barcodeFiles, tally);
        }
        // Count the number of occurrences of each unmatched sample index.
        try {
            if (NUM_THREADS == 1) {
                workers[0].call();
            } else {
                final ExecutorService threadPool = Executors.newFixedThreadPool(NUM_THREADS);
                List<Future<Boolean>> futures = threadPool.invokeAll(Arrays.stream(workers).collect(Collectors.toList()));
                for (final Future<Boolean> f: futures) {
                    if (!f.get()) {
                        throw new RuntimeException("Unpossible!");
                    }
                }
            }
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
        final long totalReads = Arrays.stream(workers).mapToLong(w -> w.totalReads).sum();
        final long totalUnmatchedReads = Arrays.stream(workers).mapToLong(w -> w.totalUnmatchedReads).sum();

        // Find the MAX_OUTPUT most frequent sample indices
        final TreeSet<IndexAndCount> indexAndCounts = new TreeSet<>();

        // Note that there for sample indices with counts that are equal, at the tail end of MAX_OUTPUT, it is
        // arbitrary which will be included in the output.
        final Collection<Map.Entry<String, LongAdder>> entries;
        if (MIN_COUNT == null) {
            entries = tally.entrySet();
        } else {
            entries = tally.entrySet().stream().filter(entry -> entry.getValue().intValue() >= MIN_COUNT).collect(Collectors.toList());
        }
        for (Map.Entry<String, LongAdder> entry: entries) {
            final int count = entry.getValue().intValue();
            if (MAX_OUTPUT == null || indexAndCounts.size() < MAX_OUTPUT) {
                indexAndCounts.add(new IndexAndCount(entry.getKey(), count));
            } else if (count > indexAndCounts.first().count) {
                indexAndCounts.remove(indexAndCounts.first());
                indexAndCounts.add(new IndexAndCount(entry.getKey(), count));
            }
        }

        // Create metrics beans for the most frequent sample indices, and write to file.
        MetricsFile<UnmatchedSampleIndexMetrics, Integer> outFile = getMetricsFile();
        outFile.addHeader(MetricsUtil.PCT_COMMENT);
        for (final IndexAndCount indexAndCount: indexAndCounts.descendingSet()) {
            outFile.addMetric(new UnmatchedSampleIndexMetrics(indexAndCount.index, indexAndCount.count,
                    indexAndCount.count/(double)totalUnmatchedReads,
                    indexAndCount.count/(double)totalReads));
        }

        outFile.write(OUTPUT);

        return 0;
    }

    private class Worker
            implements Callable<Boolean> {
        private final Queue<File> barcodeFiles;
        private final ConcurrentHashMap<String, LongAdder> map;
        long totalReads = 0;
        long totalUnmatchedReads = 0;

        public Worker(Queue<File> barcodeFiles, ConcurrentHashMap<String, LongAdder> map) {
            this.barcodeFiles = barcodeFiles;
            this.map = map;
        }

        @Override
        public Boolean call() throws Exception {
            File barcodeFile = null;
            while ((barcodeFile = barcodeFiles.poll()) != null) {
                LOG.info("Processing", barcodeFile);
                final TabbedInputParser parser = new TabbedInputParser(false, barcodeFile);
                for (final String[] row : parser) {
                    ++totalReads;
                    if ("N".equals(row[1])) {
                        ++totalUnmatchedReads;
                        map.computeIfAbsent(row[0], k -> new LongAdder()).increment();
                    }
                }
                CloserUtil.close(parser);
            }
            return true;
        }
    }

    private static class IndexAndCount
            implements Comparable<IndexAndCount> {
        final String index;
        final long count;

        public IndexAndCount(String index, long count) {
            this.index = index;
            this.count = count;
        }

        @Override
        public int compareTo(IndexAndCount o) {
            return Long.compare(this.count, o.count);
        }
    }

    private static class BarcodeFileFilter
            implements FileFilter {
        @Override
        public boolean accept(File pathname) {
            if (pathname.isFile()) {
                final String name = pathname.getName();
                return (name.endsWith("_barcode.txt.gz") || name.endsWith("_barcode.txt"));
            } else {
                return false;
            }
        }
    }
}
