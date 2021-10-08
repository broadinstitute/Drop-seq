package org.broadinstitute.dropseqrna.sbarro;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.Sbarro;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.PrintStream;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.StreamSupport;

@CommandLineProgramProperties(
        summary = "Creates a matrix of quality scores from a fastq.",
        oneLineSummary = "Creates a matrix of quality scores from a fastq.",
        programGroup = Sbarro.class)
public class QScoreCount extends CommandLineProgram {
    private static final Log log = Log.getInstance(QScoreCount.class);
    private static final int MAX_QSCORE = 100;
    private static final int MAX_READ_LENGTH = 450;

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "A fastq file.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Where to write the results.")
    public File OUTPUT;

    @Argument(shortName = "S", doc = "How many reads to skip.", optional = true)
    public long SKIP_READS = 0L;

    @Argument(shortName = "C", doc = "How many reads to analyze.", optional = true)
    public long COUNT_READS = -1L;

    @Argument(shortName = "P", doc = "Print progress after this many reads.", optional = true)
    public long PRINT_READ_PROGRESS = 10000;

    @Argument(shortName = "T", doc = "Number of parallel threads to use for processing.", optional = true)
    public int THREADS_PARALLEL = 1;

    @Override
    protected int doWork() {
        // Automatically debug small numbers of reads.
        if (0 <= COUNT_READS && COUNT_READS <= 20) {
            Log.setGlobalLogLevel(Log.LogLevel.DEBUG);
        }

        IOUtil.assertFileIsReadable(this.INPUT);
        IOUtil.assertFileIsWritable(this.OUTPUT);

        final IndexedIterator<FastqRecord, FastqReader> reader =
                new IndexedIterator<>(new FastqReader(INPUT), SKIP_READS, COUNT_READS);

        final ForkJoinPool forkJoinPool = new ForkJoinPool(THREADS_PARALLEL);
        final MatrixAggregator aggregator;
        try {
            aggregator = forkJoinPool.submit(() ->
                    StreamSupport
                            .stream(reader.spliterator(), true)
                            .collect(
                                    () -> new MatrixAggregator(log, PRINT_READ_PROGRESS),
                                    MatrixAggregator::accept,
                                    MatrixAggregator::combine
                            )
            ).get();
        } catch (InterruptedException | ExecutionException exception) {
            throw new RuntimeException(exception);
        }

        final PrintStream outStream = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(OUTPUT));
        for (int row = 0; row <= MAX_READ_LENGTH; row++) {
            for (int col = 0; col <= MAX_QSCORE; col++) {
                outStream.print(aggregator.countMatrix[row][col]);
                if (col != MAX_QSCORE) {
                    outStream.print('\t');
                } else {
                    outStream.println();
                }
            }
        }

        CloserUtil.close(reader);
        CloserUtil.close(outStream);

        return 0;
    }

    /**
     * Stock main method.
     */
    public static void main(final String[] args) {
        System.exit(new QScoreCount().instanceMain(args));
    }

    private static class MatrixAggregator {
        private final Log log;
        private final long printReadProgress;
        private final long[][] countMatrix = new long[MAX_READ_LENGTH + 1][MAX_QSCORE + 1];

        private MatrixAggregator(final Log log, final long printReadProgress) {
            this.log = log;
            this.printReadProgress = printReadProgress;
        }

        public void accept(final IndexedItem<FastqRecord> numberedFastqRecord) {
            final byte[] scores = numberedFastqRecord.getItem().getBaseQualities();
            for (int i = 0; i < scores.length; i++) {
                int score = scores[i];
                this.countMatrix[Math.min(i, MAX_READ_LENGTH)][Math.min(score, MAX_QSCORE)]++;
            }

            if (Log.isEnabled(Log.LogLevel.DEBUG)) {
                log.debug("read " + numberedFastqRecord.getIndex());
                log.debug(Arrays.toString(numberedFastqRecord.getItem().getBaseQualities()));
                log.debug();
            }

            if (numberedFastqRecord.getIndex() % printReadProgress == 0 && numberedFastqRecord.getIndex() > 0) {
                log.info(String.format("Processed read %,d", numberedFastqRecord.getIndex()));
            }
        }

        public void combine(final MatrixAggregator other) {
            for (int row = 0; row <= MAX_READ_LENGTH; row++) {
                for (int col = 0; col <= MAX_QSCORE; col++) {
                    this.countMatrix[row][col] += other.countMatrix[row][col];
                }
            }
        }
    }
}
