package org.broadinstitute.dropseqrna.sbarro;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
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
import org.broadinstitute.dropseqrna.utils.PredicateFilteredIterator;
import org.broadinstitute.dropseqrna.utils.ProgressLoggingIterator;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.TabbedInputParser;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

@CommandLineProgramProperties(
        summary = "Output reads that match a CBC and UMI that previously detected repeats, but otherwise do not match.",
        oneLineSummary = "Output reads that match a CBC and UMI that previously detected repeats, but otherwise do not match.",
        programGroup = Sbarro.class)
public class HttCbcSelector extends CommandLineProgram {
    private final Log log = Log.getInstance(HttCbcSelector.class);

    @Argument(
            shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "A bam file with cellular and molecular barcode tags."
    )
    public File INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output bam path.")
    public File OUTPUT;

    @Argument(shortName = "A", doc = "Allow tsv list of allowed CBC/UMI.")
    public File ALLOWED_CBC;

    @Argument(shortName = "F", doc = "Filter reads but do not add tags", optional = true)
    public boolean FILTER_ONLY = false;

    @Argument(shortName = "S", doc = "How many reads to skip.", optional = true)
    public long SKIP_READS = 0L;

    @Argument(shortName = "C", doc = "How many reads to analyze.", optional = true)
    public long COUNT_READS = -1L;

    @Argument(shortName = "P", doc = "Print progress after this many reads.", optional = true)
    public int PRINT_READ_PROGRESS = 100000;

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
        IOUtil.assertFileIsWritable(this.OUTPUT);

        log.info("Reading existing cell barcode umis");
        final TabbedInputParser strings;
        try {
            strings = new TabbedInputParser(false, this.ALLOWED_CBC.getCanonicalFile());
        } catch (IOException ioException) {
            throw new RuntimeException(ioException);
        }
        final Map<String, Set<String>> cellBarcodeUmis = strings
                .stream()
                .collect(
                        Collectors.groupingBy(
                                x -> x[0],
                                Collectors.mapping(
                                        x -> x[1],
                                        Collectors.toSet()
                                )
                        )
                );
        log.info("Done reading existing cell barcode umis");

        final SamReader inputSam = SamReaderFactory.makeDefault().open(INPUT);
        final SAMFileWriter outputBam =
                new SAMFileWriterFactory()
                        .setCreateIndex(false)
                        .makeSAMOrBAMWriter(inputSam.getFileHeader(), true, this.OUTPUT);

        final Iterator<SAMRecord> iterator =
                new PredicateFilteredIterator<>(
                        new ProgressLoggingIterator(
                                inputSam.iterator(),
                                new ProgressLogger(this.log, PRINT_READ_PROGRESS)
                        ),
                        new ExistingCbcUmiPredicate(cellBarcodeUmis)
                );

        final IndexedIterator<SAMRecord, Iterator<SAMRecord>> reader =
                new IndexedIterator<>(iterator, SKIP_READS, COUNT_READS);

        if (FILTER_ONLY) {
            StreamSupport
                    .stream(reader.spliterator(), false)
                    .map(IndexedItem::getItem)
                    .forEach(outputBam::addAlignment);
        } else {
            StreamSupport
                    .stream(reader.spliterator(), false)
                    .map(this::processRecord)
                    .forEach(x ->
                            {
                                final SAMRecord read = x.numberedSAMRecord.getItem();
                                if (x.skipReason != null) {
                                    read.setAttribute("XK", x.skipReason);
                                } else {
                                    read.setAttribute("XP", x.glutamineCount);
                                }
                                outputBam.addAlignment(read);
                            }
                    );
        }

        CloserUtil.close(inputSam);
        CloserUtil.close(outputBam);

        return 0;
    }

    private HttRecordSummary processRecord(final IndexedItem<SAMRecord> indexedSAMRecord) {
        final String read = indexedSAMRecord.getItem().getReadString().toUpperCase(Locale.ROOT);

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
                    indexedSAMRecord,
                    "PREFIX",
                    -1
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
                    indexedSAMRecord,
                    "SUFFIX",
                    -1
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
                    indexedSAMRecord,
                    "TRIPLETS",
                    -1
            );
        }

        return new HttRecordSummary(
                indexedSAMRecord,
                null,
                glutamineCount
        );
    }

    /**
     * Stock main method.
     */
    public static void main(final String[] args) {
        System.exit(new HttCbcSelector().instanceMain(args));
    }

    @SuppressWarnings("ClassCanBeRecord")
    private static class HttRecordSummary {
        private final IndexedItem<SAMRecord> numberedSAMRecord;
        private final String skipReason;
        private final int glutamineCount;

        private HttRecordSummary(
                final IndexedItem<SAMRecord> numberedSAMRecord,
                final String skipReason,
                final int glutamineCount
        ) {
            this.numberedSAMRecord = numberedSAMRecord;
            this.skipReason = skipReason;
            this.glutamineCount = glutamineCount;
        }
    }
}
