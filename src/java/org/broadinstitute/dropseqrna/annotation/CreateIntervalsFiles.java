/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2015 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.dropseqrna.annotation;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import org.broadinstitute.dropseqrna.cmdline.MetaData;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.File;

@CommandLineProgramProperties(
        usage = "Create standard Drop-seq intervals files: consensus_introns, genes, rRNA, exons, intergenic",
        usageShort = "Create standard Drop-seq intervals files",
        programGroup = MetaData.class
)public class CreateIntervalsFiles extends CommandLineProgram {


    @Option(shortName = StandardOptionDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, doc="The reference sequence dictionary.")
    public File SEQUENCE_DICTIONARY;

    @Option(doc="Gene definitions used to generate intervals files")
    public File REDUCED_GTF;

    @Option(shortName=StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Directory where intervals files are written")
    public File OUTPUT;

    @Option(doc="intervals files are named using this prefix")
    public String PREFIX;

    public static void main(final String[] args) {
        new CreateIntervalsFiles().instanceMainWithExit(args);
    }

    private enum ReducedGtfColumn {
        chr, start, end, annotationType, transcriptType, strand, gene_name
    }

    private enum AnnotationType {
        gene, exon, consensus_intron
    }

    private enum TranscriptType {
        rRNA
    }

    @Override
    protected int doWork() {
        IOUtil.assertDirectoryIsWritable(OUTPUT);

        final SAMFileHeader samHeader = SamReaderFactory.makeDefault().open(SEQUENCE_DICTIONARY).getFileHeader();
        final IntervalList genes = new IntervalList(samHeader);
        final IntervalList exons = new IntervalList(samHeader);
        final IntervalList consensusIntrons = new IntervalList(samHeader);
        final IntervalList rRNA = new IntervalList(samHeader);

        final TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(REDUCED_GTF);
        for (TabbedTextFileWithHeaderParser.Row row : parser) {
            try {
                final String chr = row.getField(ReducedGtfColumn.chr.name());
                final int start = Integer.parseInt(row.getField(ReducedGtfColumn.start.name()));
                final int end = Integer.parseInt(row.getField(ReducedGtfColumn.end.name()));
                final String annotationType = row.getField(ReducedGtfColumn.annotationType.name());
                final String transcriptType = row.getField(ReducedGtfColumn.transcriptType.name());
                final String strand = row.getField(ReducedGtfColumn.strand.name());
                final String name;

                if (annotationType.equals(AnnotationType.gene.name())) {
                    name = row.getField(ReducedGtfColumn.gene_name.name());
                } else {
                    name = makeIntervalName(chr, start, end);
                }
                final Interval interval = new Interval(chr, start, end, strand.equals("-"), name);

                if (annotationType.equals(AnnotationType.gene.name())) {
                    genes.add(interval);
                    if (transcriptType.contains(TranscriptType.rRNA.name())) {
                        // e.g. transcriptType "Mt_rRNA" is considered rRNA
                        rRNA.add(interval);
                    }
                } else if (annotationType.equals(AnnotationType.exon.name())) {
                    exons.add(interval);
                } else if (annotationType.equals(AnnotationType.consensus_intron.name())) {
                    consensusIntrons.add(interval);
                }
            } catch (NumberFormatException e) {
                throw new RuntimeException(String.format("Bad numeric value in file %s, line:\n %s",
                        REDUCED_GTF.getAbsolutePath(), row.getCurrentLine()), e);
            }
        }
        CloserUtil.close(parser);

        final IntervalList genome = new IntervalList(samHeader);
        for (final SAMSequenceRecord sequenceRecord : samHeader.getSequenceDictionary().getSequences()) {
            genome.add(new Interval(sequenceRecord.getSequenceName(), 1, sequenceRecord.getSequenceLength()));
        }

        final IntervalList inverted = IntervalList.invert(genes);
        final IntervalList intergenic = new IntervalList(samHeader);
        for (final Interval interval : inverted.getIntervals()) {
            intergenic.add(new Interval(interval.getContig(), interval.getStart(), interval.getEnd(), false,
                    makeIntervalName(interval.getContig(), interval.getStart(), interval.getEnd())));
        }
        write(genes, "genes");
        write(exons, "exons");
        write(consensusIntrons, "consensus_introns");
        write(rRNA, "rRNA");
        write(intergenic, "intergenic");

        return 0;
    }

    private String makeIntervalName(final String chr, final int start, final int end) {
        return String.format("%s:%d-%d", chr, start, end);
    }

    private void write(final IntervalList list, final String intervalType) {
        list.sorted().write(makeIntervalFile(intervalType));
    }

    private File makeIntervalFile(final String intervalType) {
        return new File(OUTPUT, PREFIX + "." + intervalType + ".intervals");
    }
}
