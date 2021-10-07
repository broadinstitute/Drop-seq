package org.broadinstitute.dropseqrna.sbarro;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
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
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Locale;
import java.util.PriorityQueue;
import java.util.Queue;

@CommandLineProgramProperties(
        summary = "Work in progress to search for CAG repeats in HTT.",
        oneLineSummary = "Work in progress to search for CAG repeats in HTT.",
        programGroup = Sbarro.class)
public class HttRepeats extends CommandLineProgram {
    private final Log log = Log.getInstance(HttRepeats.class);

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "A fastq file.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Where to write the results.")
    public File OUTPUT;

    @Argument(shortName = "S", doc = "How many reads to skip.", optional = true)
    public long SKIP_READS = 0L;

    @Argument(shortName = "C", doc = "How many reads to analyze.", optional = true)
    public long COUNT_READS = Long.MAX_VALUE;

    @Argument(shortName = "Q", doc = "How many Q2 tail lengths to output.", optional = true)
    public long PRINT_READ_Q2S = 10;

    @Argument(shortName = "P", doc = "Print progress after this many reads.", optional = true)
    public long PRINT_READ_PROGRESS = 1000;

    @SuppressWarnings("CommentedOutCode")
    @Override
    protected int doWork() {
        // Automatically debug small numbers of reads.
        if (COUNT_READS <= 20) {
            Log.setGlobalLogLevel(Log.LogLevel.DEBUG);
        }

        IOUtil.assertFileIsReadable(this.INPUT);
        IOUtil.assertFileIsWritable(this.OUTPUT);

        final FastqReader reader = new FastqReader(INPUT);

        // Sequence that should appear before and after mostly CAG repeats, with possibly a CAA interspersed.
        final String cagPrefix = "CCTTCGAGTCCCTCAAGTCCTTC";
        final String cagSuffix = "CCGCCACCG";
        final FindSubSequence prefixAligner = new FindSubSequence(cagPrefix);
        final FindSubSequence suffixAligner = new FindSubSequence(cagSuffix);

        final PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(OUTPUT));

        long readNumber = 0L;
        while (reader.hasNext() && readNumber < SKIP_READS) {
            reader.next();
            readNumber++;
        }

        // Count how many reads we can't find a good match for the cagPrefix and thus we skip the read.
        long skippedReads = 0;
        // For debugging keep track of the first skipped read number.
        long firstSkipped = -1;
        // When we search for the suffix, have we traversed a multiple of 3 bases? If not how far are we shifted, 0/1/2?
        final long[] shiftCounts = new long[3];
        // Calculate the average length of the Illumina Q2 tails.
        final Mean countQ2sAvg = new Mean();
        // Save a queue of the read indicies with the longes Q2 tails.
        final Queue<ReadQ2s> longestReadQ2s = new PriorityQueue<>(readQ2sComparator);

        out.println("READ_NUM\tTRIPLET_NUM\tBASES\tBASELINE\tCAG\tCAA\tCCG\tCCA");
        while (reader.hasNext() && readNumber < SKIP_READS + COUNT_READS) {
            if (readNumber % PRINT_READ_PROGRESS == 0 && readNumber > 0) {
                log.info(String.format("Processing read %,d", readNumber));
            }

            final FastqRecord record = reader.next();
            final String read = record.getReadString().toUpperCase(Locale.ROOT);
            final byte[] phredScores = record.getBaseQualities();
//            final String read = "CCTTCGAGTCCCTCAAGTCCTTCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAACAGCCGCCACCGCCGCCGCCGCCGCCGCCGCCTCCTCAGCTTCCTCAGCCGCCGCCGCAGGCACAGCCGCTGCTGCCTCAGCCGCAGCCGCCCCCGCCGCCGCCCCCGCCGCCACCCGGCCCGGCTGTGGCTGAGGAGCCGCTGCACCGACAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
//            final String qualities = "3>>ABAAFFFFFGFGGGGGGGGHHFG4GHHHHHHHHHHFHFHFHFHHHHHHHHHHHFFHFFHGHHHHHGGHHHHFHGGGGGHGGGGGGEGGGFGGCGGFGGGGGGEFHHHHHHHBFHHHGHCDFD?GG@DFADAFGGGGFAGFGFFFAFGGGFFFFFFF?CFF=DFFCFFFFFFFAFF;CFFFFDDFDFFCFFFFDFFFFFFFE.FFFD?ED--9.BBFFF--;9A?BFFFFACFFFFFCC;BFFFCFA=";
//            final byte[] phredScores = htsjdk.samtools.SAMUtils.fastqToPhred(qualities);

            final SubSequenceResultI prefixSubSeq = prefixAligner.findSequenceLocalAlignment(read);

            /*
            Skip the read if:
            - Too many mismatches
            - Prefix wasn't found near the start
             */
            final int prefixEdits = prefixSubSeq.getEditDistance().getEditDistance();
            final int prefixStart = prefixSubSeq.getStart();
            if (prefixEdits > 1 || prefixStart > 2) {
                if (firstSkipped < 0) {
                    firstSkipped = readNumber;
                }
                skippedReads++;
                readNumber++;
                continue;
            }

            // Find where the suffix starts, searching after the end of the prefix.
            final String readPostPrefix = read.substring(prefixSubSeq.getEnd());
            final SubSequenceResultI suffixSubSeq = suffixAligner.findSequenceLocalAlignment(readPostPrefix);
            // Calculate the length between the read between the end of the prefix and the beginning of the suffix.
            // Store only the count of frameshifts.
            final int shift = (suffixSubSeq.getStart() - 1) % 3;
            shiftCounts[shift]++;

            // Count the number of scores == Q2 on the tail of the read.
            int countQ2s = 0;
            for (int baseNumber = read.length() - 1; baseNumber >= 0; baseNumber--) {
                if (phredScores[baseNumber] != 2) {
                    break;
                }
                countQ2s++;
            }
            if (countQ2s > 0) {
                countQ2sAvg.increment(countQ2s);

                final ReadQ2s readQ2s = new ReadQ2s(readNumber, countQ2s);

                // Check if either: a) we need to store more reads or b) we need to kick out the lowest Q2.
                if (longestReadQ2s.size() < PRINT_READ_Q2S ||
                        readQ2sComparator.compare(longestReadQ2s.peek(), readQ2s) < 0) {
                    // Are we full? If so remove the read with the fewest number of Q2s.
                    if (longestReadQ2s.size() == PRINT_READ_Q2S) {
                        longestReadQ2s.remove();
                    }
                    // Add the new read to the queue.
                    longestReadQ2s.add(readQ2s);
                }
            }

            // Loop through complete triplets.
            final int tripletCount = readPostPrefix.length() / 3;
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

                out.printf(
                        "%d\t%d\t%s\t%f\t%f\t%f\t%f\t%f%n",
                        readNumber,
                        tripletNumber,
                        bases,
                        baselineScore,
                        cagScore,
                        caaScore,
                        ccgScore,
                        ccaScore
                );
            }

            if (Log.isEnabled(Log.LogLevel.DEBUG)) {
                log.debug("read " + readNumber);
                log.debug(String.format("%-50s%-50s%s", toString(suffixSubSeq), toString(prefixSubSeq), read));
                log.debug(Arrays.toString(record.getBaseQualities()));
                log.debug(read);
                log.debug();
            }

            readNumber++;
        }

        CloserUtil.close(reader);
        CloserUtil.close(out);

        log.info("Done!");
        log.info("Skipped reads: " + skippedReads);
        log.info("First skipped read (zero based): " + firstSkipped);
        log.info("Repeats shifts: " + Arrays.toString(shiftCounts));
        if (PRINT_READ_Q2S > 0) {
            log.info("Top Q2 tails (index: count):");
            // Priority queues are not sorted by default.
            longestReadQ2s.stream().sorted(readQ2sComparator.reversed()).forEach(
                    readQ2s -> log.info(String.format("  %d: %d", readQ2s.readNumber, readQ2s.countQ2s))
            );
        }
        log.info(String.format("Avg Q2 tails length: %.2f (N=%d)", countQ2sAvg.getResult(), countQ2sAvg.getN()));

        return 0;
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
    private static String toString(final SubSequenceResultI result) {
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

    /**
     * Store a read number and the count of Q2 quality scores on the read tail.
     */
    @SuppressWarnings("ClassCanBeRecord")
    private static class ReadQ2s {
        final long readNumber;
        final int countQ2s;

        private ReadQ2s(final long readNumber, final int countQ2s) {
            this.readNumber = readNumber;
            this.countQ2s = countQ2s;
        }
    }

    /**
     * Sort ReadQ2 instances by Q2 counts ascending, then ordering by read number descending.
     */
    private static final Comparator<ReadQ2s> readQ2sComparator = (o1, o2) -> {
        final int countDiff = o1.countQ2s - o2.countQ2s;
        return countDiff != 0 ? countDiff : (int) (o2.readNumber - o1.readNumber);
    };
}
