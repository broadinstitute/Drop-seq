package org.broadinstitute.dropseqrna.sbarro;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.Sbarro;
import org.broadinstitute.dropseqrna.sbarro.utils.FindSubSequence;
import org.broadinstitute.dropseqrna.sbarro.utils.SubSequenceResultI;
import org.broadinstitute.dropseqrna.utils.ProgressLoggingIterator;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.CellBarcodeFilteringIterator;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.TabbedInputParser;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

@CommandLineProgramProperties(
        summary = "Work in progress to search for CAG repeats in HTT.",
        oneLineSummary = "Work in progress to search for CAG repeats in HTT.",
        programGroup = Sbarro.class)
public class HttTripletGen extends CommandLineProgram {
    private final Log log = Log.getInstance(HttTripletGen.class);

    @Argument(
            shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "A bam file with cellular and molecular barcode tags."
    )
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Where to write the results.")
    public File OUTPUT;

    @Argument(shortName = "A", doc = "Allow gziped list of allowed CBC barcodes", optional = true)
    public File ALLOWED_CBC = null;

    @Argument(shortName = "S", doc = "How many reads to skip.", optional = true)
    public long SKIP_READS = 0L;

    @Argument(shortName = "C", doc = "How many reads to analyze.", optional = true)
    public long COUNT_READS = -1L;

    @Argument(shortName = "P", doc = "Print progress after this many reads.", optional = true)
    public int PRINT_READ_PROGRESS = 100000;

    @Argument(shortName = "T", doc = "Number of parallel threads to use for processing.", optional = true)
    public int THREADS_PARALLEL = 1;

    // Sequence that should appear before mostly CAG repeats, with possibly a CAA towards the end.
    private static final String cagPrefix = "CCTTCGAGTCCCTCAAGTCCTTC";

    private static final String CELL_BARCODE_TAG = "XC";
    private static final String MOLECULAR_BARCODE_TAG = "XM";

    @Override
    protected int doWork() {
        // Automatically debug small numbers of reads.
        if (0 <= COUNT_READS && COUNT_READS <= 20) {
            Log.setGlobalLogLevel(Log.LogLevel.DEBUG);
        }

        IOUtil.assertFileIsReadable(this.INPUT);
        IOUtil.assertFileIsWritable(this.OUTPUT);

        final SamReader inputSam = SamReaderFactory.makeDefault().open(INPUT);

        final ProgressLoggingIterator progressIterator = new ProgressLoggingIterator(
                inputSam.iterator(),
                new ProgressLogger(this.log, PRINT_READ_PROGRESS)
        );

        final Iterator<SAMRecord> cbcIterator;
        if (this.ALLOWED_CBC == null) {
            cbcIterator = progressIterator;
        } else {
            log.info("Reading allowed cell barcodes");
            final TabbedInputParser strings;
            try {
                strings = new TabbedInputParser(false, this.ALLOWED_CBC.getCanonicalFile());
            } catch (IOException ioException) {
                throw new RuntimeException(ioException);
            }
            final Set<String> cellBarcodes = strings.stream().map(tokens -> tokens[0]).collect(Collectors.toSet());
            log.info("Done reading allowed cell barcodes");
            cbcIterator = new CellBarcodeFilteringIterator(progressIterator, CELL_BARCODE_TAG, cellBarcodes);
        }

        final IndexedIterator<SAMRecord, Iterator<SAMRecord>> reader =
                new IndexedIterator<>(cbcIterator, SKIP_READS, COUNT_READS);

        final ForkJoinPool forkJoinPool = new ForkJoinPool(THREADS_PARALLEL);
        final HttRecordAggregator aggregator;
        try {
            aggregator = forkJoinPool.submit(() ->
                    StreamSupport
                            .stream(reader.spliterator(), true)
                            .map(this::processRecord)
                            .collect(
                                    HttRecordAggregator::new,
                                    HttRecordAggregator::accept,
                                    HttRecordAggregator::combine
                            )
            ).get();
        } catch (InterruptedException | ExecutionException exception) {
            throw new RuntimeException(exception);
        }

        final PrintStream glutamineRepeatsOut =
                new ErrorCheckingPrintStream(IOUtil.openFileForWriting(this.OUTPUT));
        glutamineRepeatsOut.println("CBC\tUMI\tREPEAT_LENGTH\tREAD_COUNT");
        aggregator.umiRepeatCounts.entrySet().stream().sorted(Map.Entry.comparingByKey()).forEach(cbcUmiEntry ->
                {
                    final String cbc = cbcUmiEntry.getKey().cbc;
                    final String umi = cbcUmiEntry.getKey().umi;
                    cbcUmiEntry.getValue().entrySet().stream().sorted(Map.Entry.comparingByKey()).forEach(countEntry ->
                            glutamineRepeatsOut
                                    .printf("%s\t%s\t%d\t%d%n", cbc, umi, countEntry.getKey(), countEntry.getValue())
                    );
                }
        );

        CloserUtil.close(inputSam);
        CloserUtil.close(glutamineRepeatsOut);

        log.info("Done!");
        return 0;
    }

    private HttRecordSummary processRecord(final IndexedItem<SAMRecord> indexedSAMRecord) {
        final String read = indexedSAMRecord.getItem().getReadString().toUpperCase(Locale.ROOT);

        final FindSubSequence prefixAligner = new FindSubSequence(cagPrefix);

        /*
        Skip the read if:
        - Too many mismatches
        - Prefix wasn't found near the start
         */
        final SubSequenceResultI prefixSubSeq = prefixAligner.findSequenceLocalAlignment(read);
        final int prefixEdits = prefixSubSeq.getEditDistance().getEditDistance();
        final int prefixStart = prefixSubSeq.getStart();
        if (prefixEdits > 1 || prefixStart > 2) {
            return new HttRecordSummary(
                    indexedSAMRecord,
                    0
            );
        }

        // Find where the suffix starts, searching after the end of the prefix.
        final char[] readPostPrefix = read.substring(prefixSubSeq.getEnd()).toCharArray();
        final char[] tripletArray = new char[3];
        // glutamines
        final char[] CAG_ARRAY = "CAG".toCharArray();
        final char[] CAA_ARRAY = "CAA".toCharArray();
        // prolines
        final char[] CCA_ARRAY = "CCA".toCharArray();
        final char[] CCC_ARRAY = "CCC".toCharArray();
        final char[] CCG_ARRAY = "CCG".toCharArray();
        final char[] CCT_ARRAY = "CCT".toCharArray();
        int otherCodeCount = 0;
        int prolineCount = 0;
        int glutamineCount = 0;
        for (int tripletNum = 0; tripletNum * 3 < readPostPrefix.length; tripletNum++) {
            System.arraycopy(readPostPrefix, tripletNum * 3, tripletArray, 0, 3);
            if (Arrays.equals(tripletArray, CAG_ARRAY) || Arrays.equals(tripletArray, CAA_ARRAY)) {
                glutamineCount++;
                // Since we saw a glutamine again, reset the proline count
                glutamineCount += prolineCount;
                prolineCount = 0;
            } else if (
                    // In order of expected frequency
                    Arrays.equals(tripletArray, CCG_ARRAY) ||
                            Arrays.equals(tripletArray, CCA_ARRAY) ||
                            Arrays.equals(tripletArray, CCC_ARRAY) ||
                            Arrays.equals(tripletArray, CCT_ARRAY)
            ) {
                prolineCount++;
            } else {
                otherCodeCount++;
            }
            if (prolineCount == 3 || otherCodeCount > 1) {
                break;
            }
        }

        // Skip the read if we didn't find any glutamines or if we found other more than one other code in the repeats.
        if (glutamineCount == 0 || otherCodeCount > 1) {
            return new HttRecordSummary(
                    indexedSAMRecord,
                    0
            );
        }

        return new HttRecordSummary(
                indexedSAMRecord,
                glutamineCount + otherCodeCount
        );
    }

    /**
     * Stock main method.
     */
    public static void main(final String[] args) {
        System.exit(new HttTripletGen().instanceMain(args));
    }

    @SuppressWarnings("ClassCanBeRecord")
    private static class HttRecordSummary {
        private final IndexedItem<SAMRecord> numberedSAMRecord;
        private final int repeatCount;

        private HttRecordSummary(
                final IndexedItem<SAMRecord> numberedSAMRecord,
                final int repeatCount
        ) {
            this.numberedSAMRecord = numberedSAMRecord;
            this.repeatCount = repeatCount;
        }
    }

    private static class HttRecordAggregator {
        private final Map<CbcUmi, Map<Integer, Long>> umiRepeatCounts = new HashMap<>();

        public void accept(final HttRecordSummary summary) {
            final String cbc = summary.numberedSAMRecord.getItem().getStringAttribute(CELL_BARCODE_TAG);
            final String umi = summary.numberedSAMRecord.getItem().getStringAttribute(MOLECULAR_BARCODE_TAG);
            if (summary.repeatCount > 0) {
                final Map<Integer, Long> counts = new HashMap<>();
                counts.put(summary.repeatCount, 1L);
                this.umiRepeatCounts.merge(
                        new CbcUmi(cbc, umi),
                        counts,
                        (existing, additions) -> {
                            additions.forEach((countRepeats, countReads) ->
                                    existing.merge(
                                            countRepeats,
                                            countReads,
                                            Long::sum
                                    )
                            );
                            return existing;
                        }
                );
            }
        }

        public void combine(final HttRecordAggregator other) {
            other.umiRepeatCounts.forEach((cbcumi, counts) ->
                    this.umiRepeatCounts.merge(
                            cbcumi,
                            counts,
                            (existing, additions) -> {
                                additions.forEach((countRepeats, countReads) ->
                                        existing.merge(
                                                countRepeats,
                                                countReads,
                                                Long::sum
                                        )
                                );
                                return existing;
                            }
                    )
            );
        }
    }

    @SuppressWarnings("ClassCanBeRecord")
    private static class CbcUmi implements Comparable<CbcUmi> {
        private final String cbc;
        private final String umi;

        private CbcUmi(final String cbc, final String umi) {
            this.cbc = cbc;
            this.umi = umi;
        }

        public String getCbc() {
            return cbc;
        }

        public String getUmi() {
            return umi;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            final CbcUmi other = (CbcUmi) o;
            return Objects.equals(this.cbc, other.cbc) && Objects.equals(this.umi, other.umi);
        }

        @Override
        public int hashCode() {
            return Objects.hash(cbc, umi);
        }

        @Override
        public int compareTo(final CbcUmi that) {
            return cbcUmiComparator.compare(this, that);
        }
    }

    private static final Comparator<CbcUmi> cbcUmiComparator =
            Comparator.comparing(CbcUmi::getCbc).thenComparing(CbcUmi::getUmi);
}
