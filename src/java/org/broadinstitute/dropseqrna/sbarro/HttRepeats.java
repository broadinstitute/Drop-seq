package org.broadinstitute.dropseqrna.sbarro;

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
import java.io.PrintStream;
import java.util.Arrays;
import java.util.Locale;

@CommandLineProgramProperties(
        summary = "Work in progress count of repeats using local alignment to match HTT prefix and suffix.",
        oneLineSummary = "Work in progress count of repeats using local alignment to match HTT prefix and suffix.",
        programGroup = Sbarro.class)
public class HttRepeats extends CommandLineProgram {
    private final Log log = Log.getInstance(HttRepeats.class);

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "A fastq file.")
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Where to write the results.")
    public File OUTPUT;

    @Argument(shortName = "S", doc = "How many records to skip.", optional = true)
    public long skip = 0L;

    @Argument(shortName = "C", doc = "How many records to output.", optional = true)
    public long count = Long.MAX_VALUE;

    @Override
    protected int doWork() {
        Log.setGlobalLogLevel(Log.LogLevel.DEBUG);

        IOUtil.assertFileIsReadable(this.INPUT);
        IOUtil.assertFileIsWritable(this.OUTPUT);

        final FastqReader reader = new FastqReader(INPUT);

        final String prefix = "CCTTCGAGTCCCTCAAGTCCTTC";
        final String suffix = "CCGCCACCG";
        final FindSubSequence findPrefix = new FindSubSequence(prefix);
        final FindSubSequence findSuffix = new FindSubSequence(suffix);

        final PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(OUTPUT));

        long readIndex = 0L;
        while (reader.hasNext() && readIndex < skip) {
            reader.next();
            readIndex++;
        }

        out.println("READ_NUM\tPOSITION\tBASES\tAVG\tBASELINE\tCAG\tCCG\tCCA");
        while (reader.hasNext() && readIndex < skip + count) {
            final FastqRecord record = reader.next();
            final String read = record.getReadString().toUpperCase(Locale.ROOT);
            final byte[] phredScores = record.getBaseQualities();
            final SubSequenceResultI resultPrefix = findPrefix.findSequenceLocalAlignment(read);
            final SubSequenceResultI resultSuffix = findSuffix.findSequenceLocalAlignment(read);
            final int tripletCount = (read.length() - resultPrefix.getMatchLength()) / 3;

            for (int tripletIndex = 0; tripletIndex < tripletCount; tripletIndex++) {
                final int tripletOffset = resultPrefix.getMatchLength() + (tripletIndex * 3);
                final double baselineScore = scoreTriplet(read, phredScores, tripletOffset, null);
                final double cagScore = scoreTriplet(read, phredScores, tripletOffset, "CAG");
                final double ccgScore = scoreTriplet(read, phredScores, tripletOffset, "CCG");
                final double ccaScore = scoreTriplet(read, phredScores, tripletOffset, "CCA");

                final double avgPValue = averageQual(phredScores, tripletOffset);
                final String bases = read.substring(tripletOffset, tripletOffset + 3);

                out.printf(
                        "%d\t%d\t%s\t%f\t%f\t%f\t%f\t%f%n",
                        readIndex,
                        tripletIndex,
                        bases,
                        avgPValue,
                        baselineScore,
                        cagScore,
                        ccgScore,
                        ccaScore
                );
            }

            if (Log.isEnabled(Log.LogLevel.DEBUG)) {
                log.debug("read " + readIndex);
                log.debug(String.format("%-50s%-50s%s", toString(resultSuffix), toString(resultPrefix), read));
                log.debug(Arrays.toString(record.getBaseQualities()));
                log.debug(read);
                log.debug();
            }

            readIndex++;
        }

        CloserUtil.close(reader);
        CloserUtil.close(out);

        return 0;
    }

    private static double averageQual(final byte[] phredScores, final int start) {
        double sum = 0;
        for (int i = 0; i < 3; i++) {
            final double pValue = LikelihoodUtils.getInstance().phredScoreToErrorProbability(phredScores[start + i]);
            sum += pValue;
        }
        return Math.log10(sum / 3D) * -10D;
    }

    private static double scoreTriplet(
            final String read,
            final byte[] phredScores,
            final int start,
            final String triplet
    ) {
        double score = 0;
        for (int i = 0; i < 3; i++) {
            final double pValue = LikelihoodUtils.getInstance().phredScoreToErrorProbability(phredScores[start + i]);
            final double prob;
            if (triplet == null || read.charAt(start + i) == triplet.charAt(i)) {
                prob = pValue;
            } else {
                prob = 1 - pValue;
            }
            score += Math.log10(prob);
        }
        return score * -10D;
    }

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
}
