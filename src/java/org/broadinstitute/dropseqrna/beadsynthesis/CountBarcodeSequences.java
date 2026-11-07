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
package org.broadinstitute.dropseqrna.beadsynthesis;

import com.google.common.collect.Iterators;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.BaseRange;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import org.broadinstitute.dropseqrna.utils.FilterBamByTag;
import org.broadinstitute.dropseqrna.utils.ProgressLoggingIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicReference;
import java.util.concurrent.atomic.LongAdder;

@CommandLineProgramProperties(
        summary = "Produce a histogram of the unique sequences for a given {read number, base range}.",
        oneLineSummary = "Produce a histogram of the unique sequences for a given {read number, base range}.",
        programGroup = DropSeq.class
)
public class CountBarcodeSequences
        extends CommandLineProgram {
    @Argument(doc = "One of more paired-end BAMs or bam_lists.", shortName = StandardOptionDefinitions.INPUT_SHORT_NAME)
    public List<File> INPUT;

    @Argument(doc = "From which read in the pair should the barcode sequence be extracted. Must be 1 or 2.")
    public int BARCODED_READ;

    @Argument(doc = "Base range to extract, seperated by a dash.  E.g 1-4.  Can extract multiple ranges by separating them by a colon.  For example 1-4:17-22 extracts the first 4 bases, then the 17-22 bases, and glues the sequence together into a single sequence for a tag.")
    public String BASE_RANGE;

    @Argument(doc = "Histogram of barcode sequences.", shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    @Argument(doc = "List of expected barcodes, one per line.  if this is specified, only barcodes in this list will be emitted.",
            optional = true)
    public File ALLOWED_BARCODES;

    @Argument(doc = "Add this value to the observed count of every barcode in the allow list.  The effect is to make " +
            "every barcode in the allow list observed at least once.")
    public int ALLOWLIST_PSEUDOCOUNT = 0;

    @Argument(doc = "Number of threads used to read and decode the input BAMs and to count barcodes.  With 1, the calling thread does everything, as before.  With more, the threads are split between reading and counting.  Records are handed to the counting threads in batches.")
    public int NUM_THREADS = 1;

    private static final int BATCH_SIZE = 1000;

    @Hidden
    @Argument(doc = "For tuning: the number of threads that read and decode the input BAMs, at least 1.  Must be given " +
            "together with COUNTING_THREADS, and then NUM_THREADS is ignored.  By default NUM_THREADS is split " +
            "between reading and counting.", optional = true)
    public Integer READER_THREADS;

    @Hidden
    @Argument(doc = "For tuning: the number of threads that count barcodes.  With 0, the calling thread counts.  Must " +
            "be given together with READER_THREADS.", optional = true)
    public Integer COUNTING_THREADS;

    private final Log log = Log.getInstance(CountBarcodeSequences.class);

    @Override
    protected String[] customCommandLineValidation() {
        IOUtil.assertFileIsWritable(OUTPUT);
        if (ALLOWED_BARCODES != null) {
            IOUtil.assertFileIsReadable(ALLOWED_BARCODES);
        }
        final ArrayList<String> list = new ArrayList<>(1);
        if (NUM_THREADS < 1) {
            list.add("NUM_THREADS must be at least 1.");
        }
        if ((READER_THREADS == null) != (COUNTING_THREADS == null)) {
            list.add("READER_THREADS and COUNTING_THREADS must be given together.");
        } else if (READER_THREADS != null && (READER_THREADS < 1 || COUNTING_THREADS < 0)) {
            list.add("READER_THREADS must be at least 1 and COUNTING_THREADS at least 0.");
        }
        if (BARCODED_READ != 1 && BARCODED_READ != 2) {
            list.add("Invalid BARCODED_READ value: " + BARCODED_READ);
        }
        return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);
    }

    @Override
    protected int doWork() {
        INPUT = FileListParsingUtils.expandFileList(INPUT);
        if (INPUT.isEmpty()) {
            throw new RuntimeException("At least one INPUT must be specified.");
        }
        INPUT.forEach(IOUtil::assertFileIsReadable);
        List<BaseRange> baseRanges = BaseRange.parseBaseRange(this.BASE_RANGE);
        final Set<String> allowedBarcodes = readAllowedBarcodes(BaseRange.getTotalRangeSize(baseRanges));
        final ConcurrentHashMap<String, LongAdder> counts = createCounts(allowedBarcodes);

        final ThreadSplit threads = getReaderAndCountingThreads();
        final SamHeaderAndIterator headerAndIterator = SamFileMergeUtil.mergeInputPaths(
                IOUtil.filesToPaths(INPUT), false, SamReaderFactory.makeDefault(), threads.readerThreads());
        final Iterator<SAMRecord> it = new ProgressLoggingIterator(headerAndIterator.iterator, new ProgressLogger(log, 10000000));
        final CountBarcodeSequenceMetrics metrics;
        if (threads.countingThreads() == 0) {
            metrics = countRecords(it, baseRanges, allowedBarcodes != null, counts);
        } else {
            metrics = countInParallel(it, threads.countingThreads(), baseRanges, allowedBarcodes != null, counts);
        }
        CloserUtil.close(headerAndIterator.iterator);
        writeMetrics(metrics, counts);
        return 0;
    }

    /**
     * @param readerThreads threads that read and decode the BAMs
     * @param countingThreads threads that count barcodes.  With 0, the calling thread counts.
     */
    private record ThreadSplit(int readerThreads, int countingThreads) {}

    /**
     * Split the threads between reading and decoding the BAMs, and counting barcodes.  READER_THREADS and
     * COUNTING_THREADS, if given, are used as they are and NUM_THREADS is ignored.  Otherwise, with up to 3 threads
     * all of them read and the calling thread counts, and with more there is a single counting thread and the rest
     * read.  Logs the split if there is more than 1 thread.
     * <p>
     * Measured with 8 compressed, unaligned inputs (1.83 billion records, 3 million allowed barcodes), in seconds,
     * for readers + counting threads:
     * <pre>
     *   2 readers + 0 counting: 301     1 reader  + 1 counting: 549
     *   2 readers + 2 counting: 303     3 readers + 1 counting: 210-219
     *   2 readers + 4 counting: 300     3 readers + 3 counting: 210
     *   2 readers + 6 counting: 302     3 readers + 5 counting: 209
     *   1 reader  + 3 counting: 538     4 readers + 0 counting: 269
     *   2 readers + 1 counting: 301     3 readers + 0 counting: 250
     * </pre>
     * With less than 4 threads, it's more optimal to add read threads to feed the main process.  Once there are 4 total
     * threads, the main process is rate limiting, making it better to allocate a thread to the counting process.
     */
    private ThreadSplit getReaderAndCountingThreads() {
        final int readerThreads;
        final int countingThreads;
        if (READER_THREADS != null) {
            readerThreads = READER_THREADS;
            countingThreads = COUNTING_THREADS;
        } else if (NUM_THREADS <= 3) {
            readerThreads = NUM_THREADS;
            countingThreads = 0;
        } else {
            readerThreads = NUM_THREADS - 1;
            countingThreads = 1;
        }
        if (readerThreads + countingThreads > 1) {
            log.info(String.format("Using %d reader thread(s) and %d counting thread(s)%s.", readerThreads,
                    countingThreads, countingThreads == 0 ? " (the main thread counts)" : ""));
        }
        return new ThreadSplit(readerThreads, countingThreads);
    }

    /**
     * @return the barcodes in ALLOWED_BARCODES, or null if there is no allow list.
     * @throws RuntimeException if a barcode is not as long as the BASE_RANGE.
     */
    private Set<String> readAllowedBarcodes(final int totalRangeSize) {
        if (ALLOWED_BARCODES == null) {
            return null;
        }
        final Set<String> allowedBarcodes = FilterBamByTag.readValues(ALLOWED_BARCODES);
        allowedBarcodes.forEach(allowedBarcode -> {
            if (allowedBarcode.length() != totalRangeSize) {
                throw new RuntimeException(String.format("BASE_RANGES implies length of %d, but value '%s' in ALLOWED_BARCODES file %s has length %d",
                        totalRangeSize, allowedBarcode, ALLOWED_BARCODES, allowedBarcode.length()));
            }
        });
        return allowedBarcodes;
    }

    /**
     * Count in a hash map rather than directly in Histogram, which is backed by a TreeMap: with millions of distinct
     * barcodes, a TreeMap lookup per read dominates run time.  The Histogram is built once at the end, and it sorts its
     * keys, so the output is the same.  The map is concurrent because several threads may count.
     * @param allowedBarcodes if not null, the map starts with every one of these barcodes, so that a single lookup both
     *                        checks the allow list and finds the counter.
     */
    private ConcurrentHashMap<String, LongAdder> createCounts(final Set<String> allowedBarcodes) {
        if (allowedBarcodes == null) {
            return new ConcurrentHashMap<>();
        }
        final ConcurrentHashMap<String, LongAdder> counts =
                new ConcurrentHashMap<>((int) (allowedBarcodes.size() / 0.75f) + 1);
        final long initialCount = ALLOWLIST_PSEUDOCOUNT > 0 ? 1 : 0;
        allowedBarcodes.forEach(bc -> {
            final LongAdder adder = new LongAdder();
            adder.add(initialCount);
            counts.put(bc, adder);
        });
        return counts;
    }

    private void writeMetrics(final CountBarcodeSequenceMetrics metrics,
                              final ConcurrentHashMap<String, LongAdder> counts) {
        final MetricsFile<CountBarcodeSequenceMetrics, String> metricsFile = getMetricsFile();
        final Histogram<String> histogram = new Histogram<>("SEQUENCE", "COUNT");
        // Allowed barcodes that were never observed (and have no pseudocount) are not emitted, as before.
        counts.forEach((bc, count) -> {
            if (count.sum() > 0) {
                histogram.increment(bc, count.sum());
            }
        });
        metricsFile.addMetric(metrics);
        metricsFile.addHistogram(histogram);
        final Writer metricsWriter = IOUtil.openFileForBufferedWriting(OUTPUT); // Enables gzip writing
        metricsFile.write(metricsWriter);
        try {
            metricsWriter.close();
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    /**
     * Pull the records from the iterator on this thread, and count them on several threads.  Records go to the
     * counting threads in batches through a small queue, so a record is only ever used by one thread at a time.  The
     * order of the records does not matter for counting.
     */
    private CountBarcodeSequenceMetrics countInParallel(final Iterator<SAMRecord> it, final int countingThreads,
                                                        final List<BaseRange> baseRanges, final boolean hasAllowList,
                                                        final ConcurrentHashMap<String, LongAdder> counts) {
        // 1. Set up the counting threads.  Each one takes batches from the queue until it gets the end marker, and
        // returns the metrics of the batches it counted.  The queue is small, so the reader cannot get far ahead of
        // the counters.
        final BlockingQueue<List<SAMRecord>> queue = new ArrayBlockingQueue<>(2 * countingThreads);
        // Compared by identity: tells a counting thread that there are no more batches.
        final List<SAMRecord> end = new ArrayList<>();
        // The first failure of a counting thread.  The others then skip their batches, so the queue keeps draining
        // and this thread is never left waiting on a full queue.
        final AtomicReference<RuntimeException> failure = new AtomicReference<>();
        final ExecutorService executor = Executors.newFixedThreadPool(countingThreads, r -> {
            final Thread t = new Thread(r, "CountBarcodeSequences-counter");
            t.setDaemon(true);
            return t;
        });
        final List<Future<CountBarcodeSequenceMetrics>> futures = new ArrayList<>(countingThreads);
        for (int i = 0; i < countingThreads; i++) {
            futures.add(executor.submit(() -> {
                final CountBarcodeSequenceMetrics threadMetrics = new CountBarcodeSequenceMetrics();
                for (List<SAMRecord> batch = take(queue); batch != end; batch = take(queue)) {
                    if (failure.get() != null) {
                        continue;
                    }
                    try {
                        add(threadMetrics, countRecords(batch.iterator(), baseRanges, hasAllowList, counts));
                    } catch (RuntimeException e) {
                        failure.compareAndSet(null, e);
                    }
                }
                return threadMetrics;
            }));
        }
        // 2. Feed the counting threads: this thread pulls records from the reader in batches and puts them on the
        // queue, which blocks while the queue is full.  Then it sends one end marker per counting thread.
        final Iterator<List<SAMRecord>> batches = Iterators.partition(it, BATCH_SIZE);
        while (failure.get() == null && batches.hasNext()) {
            put(queue, batches.next());
        }
        for (int i = 0; i < countingThreads; i++) {
            put(queue, end);
        }
        // 3. Wait for the counting threads to finish, and merge their metrics.  The counts themselves are already
        // merged, because all the threads add to the same map.
        final CountBarcodeSequenceMetrics metrics = new CountBarcodeSequenceMetrics();
        for (final Future<CountBarcodeSequenceMetrics> future : futures) {
            try {
                add(metrics, future.get());
            } catch (InterruptedException | ExecutionException e) {
                throw new RuntimeException(e);
            }
        }
        executor.shutdown();
        if (failure.get() != null) {
            throw failure.get();
        }
        return metrics;
    }

    private static List<SAMRecord> take(final BlockingQueue<List<SAMRecord>> queue) {
        try {
            return queue.take();
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        }
    }

    private static void put(final BlockingQueue<List<SAMRecord>> queue, final List<SAMRecord> batch) {
        try {
            queue.put(batch);
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        }
    }

    private static void add(final CountBarcodeSequenceMetrics total, final CountBarcodeSequenceMetrics other) {
        total.NUM_READ_PAIRS += other.NUM_READ_PAIRS;
        total.NUM_NOT_ALLOWED += other.NUM_NOT_ALLOWED;
    }

    /** Count the barcodes of the reads of one iterator, which is read to the end. */
    private CountBarcodeSequenceMetrics countRecords(final Iterator<SAMRecord> it, final List<BaseRange> baseRanges,
                                                     final boolean hasAllowList,
                                                     final ConcurrentHashMap<String, LongAdder> counts) {
        final CountBarcodeSequenceMetrics metrics = new CountBarcodeSequenceMetrics();
        while (it.hasNext()) {
            final SAMRecord rec = it.next();
            if (!rec.getReadPairedFlag()) {
                throw new RuntimeException("Unexpected unpaired read: " + rec);
            }
            if ((rec.getFirstOfPairFlag() && (BARCODED_READ == 2)) || (rec.getSecondOfPairFlag() && (BARCODED_READ == 1))) {
                continue;
            }
            ++metrics.NUM_READ_PAIRS;
            final String barcode = BaseRange.getSequenceForBaseRange(baseRanges, rec.getReadBases());
            LongAdder count = counts.get(barcode);
            if (count == null) {
                if (hasAllowList) {
                    ++metrics.NUM_NOT_ALLOWED;
                    continue;
                }
                count = counts.computeIfAbsent(barcode, k -> new LongAdder());
            }
            count.increment();
        }
        return metrics;
    }

    public static class CountBarcodeSequenceMetrics extends MetricBase {
        public long NUM_READ_PAIRS;
        public long NUM_NOT_ALLOWED;
    }
}
