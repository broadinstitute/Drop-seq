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

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.dropseqrna.cmdline.MetaData;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.sam.CreateSequenceDictionary;

import java.io.File;
import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.Set;

@CommandLineProgramProperties(
        usage = "Validate reference fasta and GTF for use in Drop-Seq, and display sequences that appear in one but " +
                "not the other, and display all gene_biotype values (transcript types)",
        usageShort = "Validate reference fasta and GTF for use in Drop-Seq",
        programGroup = MetaData.class
)public class ValidateReference extends CommandLineProgram {

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

        final CreateSequenceDictionary sequenceDictionaryCreator = new CreateSequenceDictionary();
        final SAMSequenceDictionary sequenceDictionary = sequenceDictionaryCreator.makeSequenceDictionary(REFERENCE);
        final GTFReader gtfReader = new GTFReader(GTF, sequenceDictionary);
        // Use
        final Set<String> sequencesInReference = new LinkedHashSet<>();
        for (final SAMSequenceRecord s : sequenceDictionary.getSequences()) {
            sequencesInReference.add(s.getSequenceName());
        }

        final Collection<GeneFromGTF> geneAnnotations = gtfReader.load().getAll();
        final Set<String> sequencesInGtf = new LinkedHashSet<>();
        final Set<String> transcriptTypes = new LinkedHashSet<>();
        for (final GeneFromGTF gene : geneAnnotations) {
            sequencesInGtf.add(gene.getContig());
            transcriptTypes.add(gene.getTranscriptType());
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
