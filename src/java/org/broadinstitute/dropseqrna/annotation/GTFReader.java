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

import java.io.File;
import java.util.HashSet;
import java.util.Set;

import org.broadinstitute.dropseqrna.utils.DropSeqSamUtil;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;
import picard.annotation.AnnotationException;

/**
 * Loads gene annotations from a GTF file into an OverlapDetector<Gene>.  Discards annotations that are not
 * internally consistent, e.g. transcripts on different chromosomes or different strands.
 * This borrows heavily from RefFlatReader.  Thanks picard!
 * The big difference in GTF vs RefFlat is that refFlat defines a transcript on a single line with all exons, while GTF gives you a 1 exon per line with the gene/transcript it's associated with.
 * This forces you to read 1 or more lines to assemble a transcript properly.  Coding start/stop information is also optional, and if it exists also exists on separate lines.
 * @see <a href="http://www.sanger.ac.uk/resources/software/gff/spec.html#t_2">http://www.sanger.ac.uk/resources/software/gff/spec.html#t_2</a>
 */
public class GTFReader {

    private static final Log LOG = Log.getInstance(GTFReader.class);


    private final File gtfFlatFile;
    private final SAMSequenceDictionary sequenceDictionary;
    // Keep track of errors already reported to reduce verbosity
    private final Set<String> skippedChromosomeTranscriptDescription = new HashSet<>();
    private final Set<String> unrecognizedSequences = new HashSet<>();

    public GTFReader(final File gtfFlatFile, final SAMSequenceDictionary sequenceDictionary) {
        this.gtfFlatFile = gtfFlatFile;
        this.sequenceDictionary = sequenceDictionary;
    }

    public GTFReader (final File gtfFlatFile, final File sequenceDictionary) {
        this(gtfFlatFile, DropSeqSamUtil.loadSequenceDictionary(sequenceDictionary));
    }



    static OverlapDetector<GeneFromGTF> load(final File refFlatFile, final SAMSequenceDictionary sequenceDictionary) {
        return new GTFReader(refFlatFile, sequenceDictionary).load();
    }

    public OverlapDetector<GeneFromGTF> load() {
        final FilteringGTFParser parser = new FilteringGTFParser(gtfFlatFile);
        final GeneFromGTFBuilder geneBuilder = new GeneFromGTFBuilder(parser);
        CloserUtil.close(parser);
        final OverlapDetector<GeneFromGTF> overlapDetector = new OverlapDetector<>(0, 0);

        int longestInterval = 0;
        int numIntervalsOver1MB = 0;

        while (geneBuilder.hasNext())
			try {
                // Can throw AnnotationException
                GeneFromGTF gene = geneBuilder.next();
                overlapDetector.addLhs(gene, gene);
                if (gene.length() > longestInterval) longestInterval = gene.length();
                if (gene.length() > 1000000) ++numIntervalsOver1MB;
            } catch (AnnotationException e) {
                LOG.info(e.getMessage() + " -- skipping");
            }
        LOG.debug("Longest gene: " + longestInterval + "; number of genes > 1MB: " + numIntervalsOver1MB);
        LOG.debug("Total number of genes loaded [" + overlapDetector.getAll().size() +"]");
        return overlapDetector;
    }

    public Set<String> getUnrecognizedSequences() {
        return unrecognizedSequences;
    }

    private boolean isSequenceRecognized(final String sequence) {
        return (sequenceDictionary.getSequence(sequence) != null);
    }



    private class FilteringGTFParser extends FilteredIterator<GTFRecord> {
        private FilteringGTFParser(final File gtf) {
            super(new GTFParser(gtf, ValidationStringency.STRICT));
        }

        @Override
        public boolean filterOut(final GTFRecord rec) {
            if (!isSequenceRecognized(rec.getChromosome())) {
                unrecognizedSequences.add(rec.getChromosome());
                final String transcriptDescription = rec.getGeneName() + ":" + rec.getTranscriptName();
                if (skippedChromosomeTranscriptDescription.add(transcriptDescription + "\t" + rec.getChromosome()))
					LOG.debug("Skipping " + transcriptDescription + " due to unrecognized sequence " + rec.getChromosome());
                return true;
            } else
				return false;
        }
    }
}
