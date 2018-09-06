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

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.reference.ReferenceSequenceFileFactory;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.MetaData;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import picard.PicardException;
import picard.annotation.Gene;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.IOException;
import java.util.*;

@CommandLineProgramProperties(
        summary = "Validate reference fasta and GTF for use in Drop-Seq, and display sequences that appear in one but " +
                "not the other, and display all gene_biotype values (transcript types)",
        oneLineSummary = "Validate reference fasta and GTF for use in Drop-Seq",
        programGroup = MetaData.class
)public class ValidateReference extends CommandLineProgram {

    @Override
    protected boolean requiresReference() {
        return true;
    }

    @Argument(doc="Gene annotation file in GTF format")
    public File GTF;

    @Argument(optional = true, shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc="Write report in json format, for unit testing only.")
    public File OUTPUT;

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

    static class Messages {
        List<String> transcriptsWithNoExons = new ArrayList<>();
        String baseErrors;
        List<String> sequencesOnlyInReference = new ArrayList<>();
        List<String> sequencesOnlyInGtf = new ArrayList<>();
        List<String> geneBiotypes = new ArrayList<>();
        double fractionOfSequencesOnlyInReference;
        double fractionOfSequencesOnlyInGtf;
        double fractionOfGenomeOfSequencesOnlyInReference;
    }

    private Messages messages = new Messages();

    @Override
    protected int doWork() {
        // LinkedHashSets used to preserve insertion order, which presumably has some intuitive meaning.

        final SAMSequenceDictionary sequenceDictionary = makeSequenceDictionary(REFERENCE_SEQUENCE);
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
            for (final Gene.Transcript transcript : gene) {
                if (transcript.exons.length == 0) {
                    messages.transcriptsWithNoExons.add(String.format("Gene %s, Transcript %s on sequence %s has no exons",
                            gene.getGeneID(), transcript.name, gene.getContig()));
                }
            }
        }
        messages.geneBiotypes.addAll(transcriptTypes);

        validateReferenceBases(REFERENCE_SEQUENCE);

        messages.sequencesOnlyInReference.addAll(subtract(sequencesInReference, sequencesInGtf));
        messages.sequencesOnlyInGtf.addAll(gtfReader.getUnrecognizedSequences());

        messages.fractionOfSequencesOnlyInReference = messages.sequencesOnlyInReference.size()/(double)sequencesInReference.size();
        long sizeOfOnlyInReference = 0;
        for (final String s : messages.sequencesOnlyInReference) {
            sizeOfOnlyInReference += sequenceDictionary.getSequence(s).getSequenceLength();
        }
        messages.fractionOfGenomeOfSequencesOnlyInReference = sizeOfOnlyInReference/(double)sequenceDictionary.getReferenceLength();
        messages.fractionOfSequencesOnlyInGtf = messages.sequencesOnlyInGtf.size()/(double)sequencesInGtf.size();

        // print all the problems.
        if (messages.baseErrors != null) {
            System.err.println(messages.baseErrors);
        }

        if (!messages.transcriptsWithNoExons.isEmpty()) {
            System.out.println(messages.transcriptsWithNoExons.size() + "  transcript(s) have no exons");
            for (int i = 0; i < Math.min(100, messages.transcriptsWithNoExons.size()); ++i) {
                System.out.println(messages.transcriptsWithNoExons.get(i));
            }
        }

        System.out.println("\nSequences only in reference FASTA:");
        logCollection(messages.sequencesOnlyInReference);

        System.out.println("\nSequences only in GTF:");
        logCollection(messages.sequencesOnlyInGtf);

        System.out.println("\ngene_biotype values:");
        logCollection(messages.geneBiotypes);



        System.out.println("\nFraction of sequences only in reference FASTA: " + messages.fractionOfSequencesOnlyInReference);
        System.out.println("\n(Sum of lengths of sequences only in reference FASTA)/(size of genome): " + messages.fractionOfGenomeOfSequencesOnlyInReference);
        System.out.println("\nFraction of sequences only in GTF: " + messages.fractionOfSequencesOnlyInGtf);

        if (OUTPUT != null) {
            final Gson gson = new GsonBuilder().setPrettyPrinting().create();
            try {
                ErrorCheckingPrintStream writer = new ErrorCheckingPrintStream(OUTPUT);
                writer.print(gson.toJson(messages));
            } catch (IOException e) {
                throw new RuntimeException("Exception writing " + OUTPUT.getAbsolutePath(), e);
            }
        }

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
                    messages.baseErrors = String.format("WARNING: AT least one invalid base '%c' (decimal %d) in reference sequence named %s",
                            StringUtil.byteToChar(base), base, sequence.getName());
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
