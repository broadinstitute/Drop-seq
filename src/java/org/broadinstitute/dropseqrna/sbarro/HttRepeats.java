package org.broadinstitute.dropseqrna.sbarro;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.LikelihoodUtils;
import org.broadinstitute.dropseqrna.cmdline.Sbarro;
import org.broadinstitute.dropseqrna.sbarro.utils.FindSubSequence;
import org.broadinstitute.dropseqrna.sbarro.utils.SubSequenceResultI;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;
import java.util.Objects;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.function.ToIntFunction;
import java.util.stream.StreamSupport;

@CommandLineProgramProperties(
        summary = "Work in progress to search for CAG repeats in HTT.",
        oneLineSummary = "Work in progress to search for CAG repeats in HTT.",
        programGroup = Sbarro.class)
public class HttRepeats extends CommandLineProgram {
    private final Log log = Log.getInstance(HttRepeats.class);

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "A fastq file.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Where to write the results.")
    public File OUTPUT = Paths.get("").toFile();

    @Argument(shortName = "S", doc = "How many reads to skip.", optional = true)
    public long SKIP_READS = 0L;

    @Argument(shortName = "C", doc = "How many reads to analyze.", optional = true)
    public long COUNT_READS = -1L;

    @Argument(shortName = "P", doc = "Print progress after this many reads.", optional = true)
    public long PRINT_READ_PROGRESS = 10000;

    @Argument(shortName = "OS", doc = "Output table with scores of triplets.", optional = true)
    public boolean OUTPUT_SCORES = false;

    @Argument(shortName = "T", doc = "Number of parallel threads to use for processing.", optional = true)
    public int THREADS_PARALLEL = 1;

    // Sequence that should appear before and after mostly CAG repeats, with possibly a CAA towards the end.
    private static final String cagPrefix = "CCTTCGAGTCCCTCAAGTCCTTC";
    private static final String cagSuffix = "CCGCCACCG";

    @Override
    protected int doWork() {
        // Automatically debug small numbers of reads.
        if (0 <= COUNT_READS && COUNT_READS <= 20) {
            Log.setGlobalLogLevel(Log.LogLevel.DEBUG);
        }

        IOUtil.assertFileIsReadable(this.INPUT);
        final File dir = this.OUTPUT.getAbsoluteFile();
        if (dir.exists()) {
            IOUtil.assertDirectoryIsWritable(dir);
        } else {
            IOUtil.assertDirectoryIsWritable(dir.getParentFile());
            try {
                Files.createDirectory(dir.toPath());
            } catch (IOException ioException) {
                throw new SAMException(ioException);
            }
        }

        final IndexedIterator<FastqRecord, FastqReader> reader =
                new IndexedIterator<>(new FastqReader(INPUT), SKIP_READS, COUNT_READS);

        final PrintStream scoresOut =
                !OUTPUT_SCORES ? null : new ErrorCheckingPrintStream(
                        IOUtil.openFileForWriting(new File(dir, "scores.tsv"))
                );

        if (scoresOut != null) {
            scoresOut.println("READ_NUM\tTRIPLET_NUM\tBASES\tBASELINE\tCAG\tCAA\tCCG\tCCA");
        }
        final ForkJoinPool forkJoinPool = new ForkJoinPool(THREADS_PARALLEL);
        final HttRecordAggregator aggregator;
        try {
            aggregator = forkJoinPool.submit(() ->
                    StreamSupport
                            .stream(reader.spliterator(), true)
                            .map(this::processRecord)
                            .collect(
                                    () -> new HttRecordAggregator(log, PRINT_READ_PROGRESS, scoresOut),
                                    HttRecordAggregator::accept,
                                    HttRecordAggregator::combine
                            )
            ).get();
        } catch (InterruptedException | ExecutionException exception) {
            throw new RuntimeException(exception);
        }

        final PrintStream glutamineRepeatsOut =
                new ErrorCheckingPrintStream(
                        IOUtil.openFileForWriting(new File(dir, "glutamine_repeats.tsv"))
                );
        glutamineRepeatsOut.println(
                "PREFIX_OFFSET\tPREFIX_MISMATCHES\tSUFFIX_MISMATCHES\tGLUTAMINE_COUNT\tQ2_START\tCOUNT\tFIRST"
        );
        aggregator.glutamineRepeatCounts
                .entrySet()
                .stream()
                .sorted(Map.Entry.comparingByKey(glutamineRepeatsComparator))
                .forEachOrdered(
                        entry ->
                                glutamineRepeatsOut.printf(
                                        "%d\t%d\t%d\t%d\t%d\t%d\t%d%n",
                                        entry.getKey().prefixOffset,
                                        entry.getKey().prefixMismatches,
                                        entry.getKey().suffixMismatches,
                                        entry.getKey().glutamineCount,
                                        entry.getKey().q2Start,
                                        entry.getValue().count,
                                        entry.getValue().firstReadNumber
                                )
                );

        CloserUtil.close(reader);
        CloserUtil.close(scoresOut);
        CloserUtil.close(glutamineRepeatsOut);

        log.info("Done!");
        log.info("Skipped reads: " + aggregator.skippedReads);
        log.info("First skipped read (zero based): " + aggregator.firstSkipped);
        log.info("Repeats shifts: " + Arrays.toString(aggregator.shiftCounts));
        return 0;
    }

    @SuppressWarnings("CommentedOutCode")
    private HttRecordSummary processRecord(final IndexedItem<FastqRecord> indexedFastqRecord) {
        final String read = indexedFastqRecord.getItem().getReadString().toUpperCase(Locale.ROOT);
        final byte[] phredScores = indexedFastqRecord.getItem().getBaseQualities();
//        final String read = "CCTTCGAGTCCCTCAAGTCCTTCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAACAGCCGCCACCGCCGCCGCCGCCGCCGCCGCCTCCTCAGCTTCCTCAGCCGCCGCCGCAGGCACAGCCGCTGCTGCCTCAGCCGCAGCCGCCCCCGCCGCCGCCCCCGCCGCCACCCGGCCCGGCTGTGGCTGAGGAGCCGCTGCACCGACAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
//        final String qualities = "3>>ABAAFFFFFGFGGGGGGGGHHFG4GHHHHHHHHHHFHFHFHFHHHHHHHHHHHFFHFFHGHHHHHGGHHHHFHGGGGGHGGGGGGEGGGFGGCGGFGGGGGGEFHHHHHHHBFHHHGHCDFD?GG@DFADAFGGGGFAGFGFFFAFGGGFFFFFFF?CFF=DFFCFFFFFFFAFF;CFFFFDDFDFFCFFFFDFFFFFFFE.FFFD?ED--9.BBFFF--;9A?BFFFFACFFFFFCC;BFFFCFA=";
//        final byte[] phredScores = htsjdk.samtools.SAMUtils.fastqToPhred(qualities);

        final FindSubSequence prefixAligner = new FindSubSequence(cagPrefix);
        final FindSubSequence suffixAligner = new FindSubSequence(cagSuffix);

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
                    indexedFastqRecord,
                    true,
                    prefixSubSeq,
                    null,
                    -1,
                    null,
                    null
            );
        }

        // Find where the suffix starts, searching after the end of the prefix.
        final String readPostPrefix = read.substring(prefixSubSeq.getEnd());
        final SubSequenceResultI suffixSubSeq = suffixAligner.findSequenceLocalAlignment(readPostPrefix);
        // Calculate the length between the read between the end of the prefix and the beginning of the suffix.
        // Store only the count of frameshifts.
        final int shift = (suffixSubSeq.getStart() - 1) % 3;

        // Skip the read if we didn't find the suffix, or if we didn't find a multiple of 3 bases.
        final int suffixEdits = suffixSubSeq.getEditDistance().getEditDistance();
        if (shift != 0 || suffixEdits > 1) {
            return new HttRecordSummary(
                    indexedFastqRecord,
                    true,
                    prefixSubSeq,
                    suffixSubSeq,
                    shift,
                    null,
                    null
            );
        }

        // Get a string containing just the supposed repeats.
        final String readRepeats = readPostPrefix.substring(0, suffixSubSeq.getStart() - 1);
        final String[] triplets = new String[readRepeats.length() / 3];
        // Copy the triplets into an array.
        for (int tripletNumber = 0; tripletNumber < triplets.length; tripletNumber++) {
            triplets[tripletNumber] = readRepeats.substring(tripletNumber * 3, (tripletNumber + 1) * 3);
        }
        // Count how many glutamines we find.
        int glutamineCount = 0;
        for (final String triplet : triplets) {
            if (triplet.equals("CAG") || triplet.equals("CAA")) {
                glutamineCount++;
            }
        }
        // If the triplets weren't mostly glutamines then skip this read.
        if (glutamineCount == 0 || glutamineCount + 1 < triplets.length) {
            return new HttRecordSummary(
                    indexedFastqRecord,
                    true,
                    prefixSubSeq,
                    suffixSubSeq,
                    shift,
                    null,
                    null
            );
        }

        // Find the location of the Q2 tail
        int q2Start;
        for (q2Start = phredScores.length; q2Start > 0; q2Start--) {
            if (phredScores[q2Start - 1] != 2) {
                break;
            }
        }

        final GlutamineRepeat glutamineRepeat =
                new GlutamineRepeat(
                        prefixSubSeq.getStart() - 1,
                        prefixSubSeq.getEditDistance().getEditDistance(),
                        suffixSubSeq.getEditDistance().getEditDistance(),
                        glutamineCount,
                        q2Start
                );

        final TripletScore[] tripletScores;
        if (!OUTPUT_SCORES) {
            tripletScores = null;
        } else {
            // Loop through complete triplets.
            final int tripletCount = readPostPrefix.length() / 3;
            tripletScores = new TripletScore[tripletCount];
            for (int tripletNumber = 0; tripletNumber < tripletCount; tripletNumber++) {
                final int startBase = prefixSubSeq.getMatchLength() + (tripletNumber * 3);
                final double baselineScore = scoreTriplet(read, phredScores, startBase, null);
                // CAG / CAA = Glutamine
                final double cagScore = scoreTriplet(read, phredScores, startBase, "CAG");
                final double caaScore = scoreTriplet(read, phredScores, startBase, "CAA");
                // CCG / CCA = Proline
                final double ccgScore = scoreTriplet(read, phredScores, startBase, "CCG");
                final double ccaScore = scoreTriplet(read, phredScores, startBase, "CCA");

                final String bases = read.substring(startBase, startBase + 3);

                tripletScores[tripletNumber] =
                        new TripletScore(
                                bases,
                                baselineScore,
                                cagScore,
                                caaScore,
                                ccgScore,
                                ccaScore
                        );
            }
        }

        return new HttRecordSummary(
                indexedFastqRecord,
                false,
                prefixSubSeq,
                suffixSubSeq,
                shift,
                glutamineRepeat,
                tripletScores
        );
    }

    /**
     * Return the combined error likelihood.
     *
     * Inputs are in phred scale, outputs are in log10.
     */
    private static double scoreTriplet(
            final String read,
            final byte[] phredScores,
            final int start,
            final String triplet
    ) {
        double errorPenalty = 0;
        for (int i = 0; i < 3; i++) {
            final byte phredScore = phredScores[start + i];
            final double errorLikelihood = LikelihoodUtils.getInstance().phredScoreToErrorProbability(phredScore);
            if (triplet == null || read.charAt(start + i) == triplet.charAt(i)) {
                errorPenalty += Math.log10(1 - errorLikelihood);
            } else {
                errorPenalty += Math.log10(errorLikelihood);
            }
        }
        return errorPenalty;
    }

    /**
     * Create a printable version of a sub sequence result.
     */
    private static String subSeqToString(final SubSequenceResultI result) {
        if (result == null) {
            return "<none>";
        }
        return String.format(
                "%5s (%3s - %3s) %3s",
                result.getEditDistance().getEditDistance(),
                result.getStart(),
                result.getEnd(),
                result.getSubSequence()
        );
    }

    /**
     * Stock main method.
     */
    public static void main(final String[] args) {
        System.exit(new HttRepeats().instanceMain(args));
    }

    @SuppressWarnings("ClassCanBeRecord")
    private static class GlutamineRepeat {
        final int prefixOffset;
        final int prefixMismatches;
        final int suffixMismatches;
        final int glutamineCount;
        final int q2Start;

        private GlutamineRepeat(
                final int prefixOffset,
                final int prefixMismatches,
                final int suffixMismatches,
                final int glutamineCount,
                final int q2Start
        ) {
            this.prefixOffset = prefixOffset;
            this.prefixMismatches = prefixMismatches;
            this.suffixMismatches = suffixMismatches;
            this.glutamineCount = glutamineCount;
            this.q2Start = q2Start;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            GlutamineRepeat that = (GlutamineRepeat) o;
            return this.prefixOffset == that.prefixOffset &&
                    this.prefixMismatches == that.prefixMismatches &&
                    this.suffixMismatches == that.suffixMismatches &&
                    this.glutamineCount == that.glutamineCount &&
                    this.q2Start == that.q2Start;
        }

        @Override
        public int hashCode() {
            return Objects.hash(prefixOffset, prefixMismatches, suffixMismatches, glutamineCount, q2Start);
        }
    }

    @SuppressWarnings("ClassCanBeRecord")
    private static class GlutamineRepeatCount {
        private final long count;
        private final long firstReadNumber;

        private GlutamineRepeatCount(final long count, final long firstReadNumber) {
            this.count = count;
            this.firstReadNumber = firstReadNumber;
        }

        public static GlutamineRepeatCount initial(final long readNumber) {
            return new GlutamineRepeatCount(1L, readNumber);
        }

        public static GlutamineRepeatCount combine(final GlutamineRepeatCount o1, final GlutamineRepeatCount o2) {
            return new GlutamineRepeatCount(
                    o1.count + o2.count,
                    Math.min(o1.firstReadNumber, o2.firstReadNumber)
            );
        }
    }

    @SuppressWarnings("ClassCanBeRecord")
    private static class HttRecordSummary {
        private final IndexedItem<FastqRecord> numberedFastqRecord;
        private final boolean skipped;
        private final SubSequenceResultI prefixSubSeq;
        private final SubSequenceResultI suffixSubSeq;
        private final int shift;
        private final GlutamineRepeat glutamineRepeat;
        private final TripletScore[] tripletScores;

        private HttRecordSummary(
                final IndexedItem<FastqRecord> numberedFastqRecord,
                final boolean skipped,
                final SubSequenceResultI prefixSubSeq,
                final SubSequenceResultI suffixSubSeq,
                final int shift,
                final GlutamineRepeat glutamineRepeat,
                final TripletScore[] tripletScores
        ) {
            this.numberedFastqRecord = numberedFastqRecord;
            this.skipped = skipped;
            this.prefixSubSeq = prefixSubSeq;
            this.suffixSubSeq = suffixSubSeq;
            this.shift = shift;
            this.glutamineRepeat = glutamineRepeat;
            this.tripletScores = tripletScores;
        }
    }

    private static class HttRecordAggregator {
        private final Log log;
        private final long printReadProgress;
        private final PrintStream scoresOut;

        // Count how many reads where we can't find a good match for the cagPrefix and thus we skip the read.
        private long skippedReads = 0;
        // For debugging keep track of the first skipped read number.
        private long firstSkipped = -1;
        // When we search for the suffix, have we traversed a multiple of 3 bases? If not how far are we shifted, 0/1/2?
        private final long[] shiftCounts = new long[3];

        private final Map<GlutamineRepeat, GlutamineRepeatCount> glutamineRepeatCounts =
                new HashMap<>();

        private HttRecordAggregator(final Log log, final long printReadProgress, final PrintStream scoresOut) {
            this.log = log;
            this.printReadProgress = printReadProgress;
            this.scoresOut = scoresOut;
        }

        public void accept(final HttRecordSummary summary) {
            if (summary.skipped) {
                skippedReads++;
                if (firstSkipped < 0) {
                    firstSkipped = summary.numberedFastqRecord.getIndex();
                }
            }

            if (summary.shift >= 0) {
                shiftCounts[summary.shift]++;
            }

            if (scoresOut != null && summary.tripletScores != null) {
                for (int tripletNumber = 0; tripletNumber < summary.tripletScores.length; tripletNumber++) {
                    final TripletScore triplet = summary.tripletScores[tripletNumber];
                    scoresOut.printf(
                            "%d\t%d\t%s\t%f\t%f\t%f\t%f\t%f%n",
                            summary.numberedFastqRecord.getIndex(),
                            tripletNumber,
                            triplet.bases,
                            triplet.baselineScore,
                            triplet.cagScore,
                            triplet.caaScore,
                            triplet.ccgScore,
                            triplet.ccaScore
                    );
                }
            }

            if (summary.glutamineRepeat != null) {
                glutamineRepeatCounts.merge(
                        summary.glutamineRepeat,
                        GlutamineRepeatCount.initial(summary.numberedFastqRecord.getIndex()),
                        GlutamineRepeatCount::combine
                );
            }

            if (Log.isEnabled(Log.LogLevel.DEBUG)) {
                log.debug("read " + summary.numberedFastqRecord.getIndex());
                log.debug(String.format(
                        "%-50s%-50s",
                        subSeqToString(summary.suffixSubSeq),
                        subSeqToString(summary.prefixSubSeq)
                ));
                log.debug(Arrays.toString(summary.numberedFastqRecord.getItem().getBaseQualities()));
                log.debug(summary.numberedFastqRecord.getItem().getReadString());
                log.debug();
            }

            if (summary.numberedFastqRecord.getIndex() % printReadProgress == 0 &&
                    summary.numberedFastqRecord.getIndex() > 0) {
                log.info(String.format("Processed read %,d", summary.numberedFastqRecord.getIndex()));
            }
        }

        public void combine(final HttRecordAggregator other) {
            this.skippedReads += other.skippedReads;
            if (this.firstSkipped < 0 || this.firstSkipped > other.firstSkipped && other.firstSkipped > 0) {
                this.firstSkipped = other.firstSkipped;
            }
            for (int i = 0; i < other.shiftCounts.length; i++) {
                this.shiftCounts[i] += other.shiftCounts[i];
            }
            other.glutamineRepeatCounts.forEach((glutamineRepeat, value) ->
                    this.glutamineRepeatCounts.merge(glutamineRepeat, value, GlutamineRepeatCount::combine)
            );
        }
    }

    @SuppressWarnings("ClassCanBeRecord")
    private static class TripletScore {
        private final String bases;
        private final double baselineScore;
        private final double cagScore;
        private final double caaScore;
        private final double ccgScore;
        private final double ccaScore;

        private TripletScore(
                final String bases,
                final double baselineScore,
                final double cagScore,
                final double caaScore,
                final double ccgScore,
                final double ccaScore
        ) {
            this.bases = bases;
            this.baselineScore = baselineScore;
            this.cagScore = cagScore;
            this.caaScore = caaScore;
            this.ccgScore = ccgScore;
            this.ccaScore = ccaScore;
        }
    }

    private static final Comparator<GlutamineRepeat> glutamineRepeatsComparator =
            Comparator
                    .comparingInt((ToIntFunction<GlutamineRepeat>) glutamineRepeats -> glutamineRepeats.prefixOffset)
                    .thenComparingInt(glutamineRepeats -> glutamineRepeats.prefixMismatches)
                    .thenComparingInt(glutamineRepeats -> glutamineRepeats.suffixMismatches)
                    .thenComparingInt(glutamineRepeats -> glutamineRepeats.glutamineCount)
                    .thenComparingInt(glutamineRepeats -> glutamineRepeats.q2Start);
}
