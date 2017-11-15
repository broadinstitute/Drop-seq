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

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.dropseqrna.cmdline.MetaData;
import picard.PicardException;
import picard.annotation.Gene;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.sam.CreateSequenceDictionary;

import java.io.File;
import java.util.*;

@CommandLineProgramProperties(
        usage = "Validate reference fasta and GTF for use in Drop-Seq, and display sequences that appear in one but " +
                "not the other, and display all gene_biotype values (transcript types)",
        usageShort = "Validate reference fasta and GTF for use in Drop-Seq",
        programGroup = MetaData.class
)public class ValidateReference extends CommandLineProgram {

    @SuppressWarnings("WeakerAccess")
    @Option(shortName = StandardOptionDefinitions.REFERENCE_SHORT_NAME, doc="The reference fasta")
    public File REFERENCE;

    @Option(doc="Gene annotation file in GTF format")
    public File GTF;

    public static void main(final String[] args) {
        new ValidateReference().instanceMainWithExit(args);
    }

    // Create lookup table to determine if base in fasta is valid.
    // From http://www.bioinformatics.org/sms/iupac.html
    private static final String IUPAC_CODES = "ACGTURYSWKMBDHVN";
    private static final boolean[] IUPAC_TABLE = new boolean[256];

    static {
        for (final char c : IUPAC_CODES.toCharArray()) {
            IUPAC_TABLE[c] = true;
            IUPAC_TABLE[Character.toLowerCase(c)] = true;
        }
    }

    @Override
    protected int doWork() {
        // LinkedHashSets used to preserve insertion order, which presumably has some intuitive meaning.

        final SAMSequenceDictionary sequenceDictionary = makeSequenceDictionary(REFERENCE);
        final GTFReader gtfReader = new GTFReader(GTF, sequenceDictionary);
        // Use
        final Set<String> sequencesInReference = new LinkedHashSet<>();
        for (final SAMSequenceRecord s : sequenceDictionary.getSequences()) {
            sequencesInReference.add(s.getSequenceName());
        }

        final Collection<GeneFromGTF> geneAnnotations = gtfReader.load().getAll();
        final Set<String> sequencesInGtf = new LinkedHashSet<>();
        final Set<String> transcriptTypes = new LinkedHashSet<>();

        final List<String> transcriptsWithNoExons = new ArrayList<>();

        for (final GeneFromGTF gene : geneAnnotations) {
            sequencesInGtf.add(gene.getContig());
            transcriptTypes.add(gene.getTranscriptType());
            for (final Gene.Transcript transcript : gene) {
                if (transcript.exons.length == 0) {
                    transcriptsWithNoExons.add(String.format("Gene %s, Transcript %s on sequence %s has no exons",
                            gene.getGeneID(), transcript.name, gene.getContig()));
                }
            }
        }

        if (!transcriptsWithNoExons.isEmpty()) {
            System.out.println(transcriptsWithNoExons.size() + "  transcript(s) have no exons");
            for (int i = 0; i < Math.min(100, transcriptsWithNoExons.size()); ++i) {
                System.out.println(transcriptsWithNoExons.get(i));
            }
        }

        validateReferenceBases(REFERENCE);

        final Set<String> onlyInReference = subtract(sequencesInReference, sequencesInGtf);
        final Set<String> onlyInGtf = gtfReader.getUnrecognizedSequences();

        System.out.println("\nSequences only in reference FASTA:");
        logCollection(onlyInReference);

        System.out.println("\nSequences only in GTF:");
        logCollection(onlyInGtf);

        System.out.println("\ngene_biotype values:");
        logCollection(transcriptTypes);

        final double fractionOfSequencesOnlyInReference = onlyInReference.size()/(double)sequencesInReference.size();
        long sizeOfOnlyInReference = 0;
        for (final String s : onlyInReference) {
            sizeOfOnlyInReference += sequenceDictionary.getSequence(s).getSequenceLength();
        }
        final double fractionOfGenomeOfSequencesOnlyInReference = sizeOfOnlyInReference/(double)sequenceDictionary.getReferenceLength();
        final double fractionOfSequencesOnlyInGtf = onlyInGtf.size()/(double)sequencesInGtf.size();
        System.out.println("\nFraction of sequences only in reference FASTA: " + fractionOfSequencesOnlyInReference);
        System.out.println("\n(Sum of lengths of sequences only in reference FASTA)/(size of genome): " + fractionOfGenomeOfSequencesOnlyInReference);
        System.out.println("\nFraction of sequences only in GTF: " + fractionOfSequencesOnlyInGtf);

        return 0;
    }

    private SAMSequenceDictionary makeSequenceDictionary(final File referenceFile) {
        final ReferenceSequenceFile refSeqFile =
                ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceFile, true);
        ReferenceSequence refSeq;
        final List<SAMSequenceRecord> ret = new ArrayList<>();
        final Set<String> sequenceNames = new HashSet<>();
        while ((refSeq = refSeqFile.nextSequence()) != null) {
            if (sequenceNames.contains(refSeq.getName())) {
                throw new PicardException("Sequence name appears more than once in reference: " + refSeq.getName());
            }
            sequenceNames.add(refSeq.getName());
            ret.add(new SAMSequenceRecord(refSeq.getName(), refSeq.length()));
        }
        return new SAMSequenceDictionary(ret);
    }


    private void validateReferenceBases(File referenceFile) {
        final ReferenceSequenceFile refSeqFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceFile, true);
        ReferenceSequence sequence;
        while ((sequence = refSeqFile.nextSequence()) != null) {
            for (final byte base: sequence.getBases()) {
                if (!IUPAC_TABLE[base]) {
                    System.err.println(String.format("WARNING: AT least one invalid base '%c' (decimal %d) in reference sequence named %s",
                            StringUtil.byteToChar(base), base, sequence.getName()));
                    break;
                }
            }
        }
    }

    private static <T> Set<T> subtract(final Set<T> setToSubtractFrom, final Set<T> setToSubtract) {
        final Set<T> ret = new LinkedHashSet<>(setToSubtractFrom);
        ret.removeAll(setToSubtract);
        return ret;
    }

    private void logCollection(final Collection<String> collection) {
        if (collection.isEmpty()) {
            System.out.println("(none)");
        } else {
            for (final String s : collection) {
                System.out.println(s);
            }
        }
    }
}
