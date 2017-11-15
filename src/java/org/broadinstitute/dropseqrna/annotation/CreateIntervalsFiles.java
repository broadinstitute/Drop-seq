/*
 * MIT License
 *
 * Copyright 2017 Broad Institute
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
import java.util.ArrayList;
import java.util.List;

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

    @Option(doc="Name(s) of MT reference sequence, for creating mt.intervals file")
    public List<String> MT_SEQUENCE;

    @Option(doc="Name(s) of non-autosome reference sequences, for creating non_autosomes.intervals file")
    public List<String> NON_AUTOSOME_SEQUENCE;

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
        final IntervalList mtIntervals;
        final IntervalList nonAutosomeIntervals;

        if (MT_SEQUENCE.isEmpty()) {
            mtIntervals = null;
        } else {
            mtIntervals = new IntervalList(createSubsetSamHeader(samHeader, MT_SEQUENCE));
        }

        if (NON_AUTOSOME_SEQUENCE.isEmpty()) {
            nonAutosomeIntervals = null;
        } else {
            nonAutosomeIntervals = new IntervalList(createSubsetSamHeader(samHeader, NON_AUTOSOME_SEQUENCE));
        }


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
                    if (MT_SEQUENCE.contains(chr)) {
                        mtIntervals.add(interval);
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
        if (mtIntervals != null) {
            write(mtIntervals, "mt");
        }

        if (nonAutosomeIntervals != null) {
            write(nonAutosomeIntervals, "non_autosomes");
        }

        return 0;
    }

    private SAMFileHeader createMtSamHeader(final SAMFileHeader samHeader) {
        final SAMFileHeader mtHeader = samHeader.clone();
        final List<SAMSequenceRecord> mtSequences = new ArrayList<>(MT_SEQUENCE.size());
        for (final String mtSequence: MT_SEQUENCE) {
            final SAMSequenceRecord sequenceRecord = mtHeader.getSequence(mtSequence);
            if (sequenceRecord == null) {
                throw new RuntimeException("MT sequence '" + mtSequence + "' not found in sequence dictionary");
            }
            mtSequences.add(sequenceRecord);
        }
        mtHeader.getSequenceDictionary().setSequences(mtSequences);
        return mtHeader;
    }

    /**
     * Create a SAM header with a subset of the sequences
     * @param samHeader header to be subsetted
     * @param sequenceNames sequences to include
     * @return clone of the original header, but with only the named sequences
     */
    private static SAMFileHeader createSubsetSamHeader(final SAMFileHeader samHeader, final List<String> sequenceNames) {
        final SAMFileHeader ret = samHeader.clone();
        final List<SAMSequenceRecord> sequenceRecords = new ArrayList<>(sequenceNames.size());
        for (final String sequence: sequenceNames) {
            final SAMSequenceRecord sequenceRecord = ret.getSequence(sequence);
            if (sequenceRecord == null) {
                throw new RuntimeException("Sequence '" + sequence + "' specified on command line but not found in sequence dictionary");
            }
            sequenceRecords.add(sequenceRecord);
        }
        ret.getSequenceDictionary().setSequences(sequenceRecords);
        return ret;
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
