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

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.dropseqrna.utils.FilteredIterator;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.Log;
import picard.annotation.AnnotationException;

/**
 * Convert a collection of GTFRecords into a collection of GeneFromGTF objects.
 * If multiple gene_versions of a gene, only the latest is selected.
 * Extensive cross-check of all the records that the gene comprises.
 * NOTE: next() can throw AnnotationException
 */
public class GeneFromGTFBuilder implements Iterator<GeneFromGTF> {

    private final Log LOG = Log.getInstance(GeneFromGTFBuilder.class);

    final Iterator<Collection<GTFRecord>> gtfRecordsByGeneIterator;

    public GeneFromGTFBuilder(final Iterator<GTFRecord> gtfRecords) {
        final Map<String, Collection<GTFRecord>> gatheredByGene = gatherByGeneName(gtfRecords);
        gtfRecordsByGeneIterator = gatheredByGene.values().iterator();
    }

    @Override
    public boolean hasNext() {
        return gtfRecordsByGeneIterator.hasNext();
    }

    /**
     * NOTE: May throw AnnotationException on invalid input
     */
    @Override
    public GeneFromGTF next() {
        final Collection<GTFRecord> recordsForGene = gtfRecordsByGeneIterator.next();
        return makeGeneFromMultiVersionGTFRecords(recordsForGene);
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    private GeneFromGTF makeGeneFromMultiVersionGTFRecords(final Collection<GTFRecord> gtfRecords) {
        final Map<Integer, Collection<GTFRecord>> gtfRecordsByGeneVersion = gatherByGeneVersion(gtfRecords);
        final List<Integer> versions = new ArrayList<>(gtfRecordsByGeneVersion.keySet());
        Collections.sort(versions);
        final int highestVersion = versions.get(versions.size() - 1);
        return makeGeneWithTranscriptsFromGTFRecords(gtfRecordsByGeneVersion.get(highestVersion));
    }

    private GeneFromGTF makeGeneWithTranscriptsFromGTFRecords(final Collection<GTFRecord> gtfRecords) {
        final GeneFromGTF gene = makeGeneFromGTFRecords(gtfRecords);
        // Remove featureType==gene before making transcripts
        final Collection<GTFRecord> nonGeneGTFRecords = CollectionUtil.makeCollection(new GeneAnnotationFilter(gtfRecords.iterator()));

        final Map<String, Collection<GTFRecord>> gtfLinesByTranscript = gatherByTranscriptId(nonGeneGTFRecords);
        for (final Map.Entry<String, Collection<GTFRecord>> entry : gtfLinesByTranscript.entrySet()) {
            if (entry.getKey() == null)
				// Skip gene entries
                continue;
            addTranscriptToGeneFromGTFRecords(gene, entry.getValue());
        }

        if (!gene.iterator().hasNext())
			throw new AnnotationException("No transcript in GTF for gene " + gene.getName());

        return gene;
    }

    private GeneFromGTF makeGeneFromGTFRecords(final Collection<GTFRecord> gtfRecords) {
        GTFRecord lineOne=gtfRecords.iterator().next();

        String geneName=lineOne.getGeneName();

        final boolean transcriptNegStrand = lineOne.isNegativeStrand();

        // Figure out the extend of the gene
        int start = Integer.MAX_VALUE;
        int end = Integer.MIN_VALUE;
        final Set<String> geneIds = new HashSet<>();
        final Set<String> chromosomes = new HashSet<>();
        for (final GTFRecord r: gtfRecords) {
            start = Math.min(start, r.getStart());
            end   = Math.max(end,   r.getEnd());
            geneIds.add(r.getGeneID());
            chromosomes.add(r.getChromosome());
        }
        if (chromosomes.size() > 1)
			throw new AnnotationException("Chromosome disagreement(" + CollectionUtil.join(chromosomes, ", ") +
                    ") in GTF file for gene " + geneName);
        final GeneFromGTF gene = new GeneFromGTF(lineOne.getChromosome(), start, end, transcriptNegStrand, geneName, lineOne.getFeatureType(),
                lineOne.getGeneID(), lineOne.getTranscriptType(), lineOne.getGeneVersion());

        for (final GTFRecord gtfRecord : gtfRecords)
			validateGTFRecord(gtfRecord, gene);

        if (geneIds.size() > 1)
			throw new AnnotationException(String.format("Multiple gene IDs for gene %s: %s", geneName, CollectionUtil.join(geneIds, ", ")));

        return gene;
    }

    /**
     * Conversion from 0-based half-open to 1-based inclusive intervals is done here.
     */
    private GeneFromGTF.TranscriptFromGTF addTranscriptToGeneFromGTFRecords(final GeneFromGTF gene,
                                                                            final Collection<GTFRecord> transcriptLines) {

        final String geneName = gene.getName();

        GTFRecord lineOne=transcriptLines.iterator().next();
        final String transcriptName = lineOne.getTranscriptName();
        final String transcriptID = lineOne.getTranscriptID();
        final String transcriptDescription = geneName + ":" + transcriptName;
        String transcriptType = lineOne.getTranscriptType();

        List<Exon> exons = new ArrayList<>();

        int transcriptionStart = Integer.MAX_VALUE;
        int transcriptionEnd = Integer.MIN_VALUE;
        int codingStart = Integer.MAX_VALUE;
        int codingEnd = Integer.MIN_VALUE;

        for (GTFRecord r: transcriptLines) {
            String featureType=r.getFeatureType();
            int start = r.getStart();
            int end = r.getEnd();

            if (featureType.equals("exon")) {
                Exon e = new Exon(start,end);
                exons.add(e);
                transcriptionStart = Math.min(transcriptionStart, start);
                transcriptionEnd = Math.max(transcriptionEnd, end);

            }
            if (featureType.equals("CDS")) {
                codingStart=Math.min(codingStart, start);
                codingEnd=Math.max(codingEnd, end);
            }
        }
        Collections.sort(exons);
        if (codingStart==Integer.MAX_VALUE) codingStart=transcriptionStart;
        if (codingEnd==Integer.MIN_VALUE) codingEnd=transcriptionEnd;

        final GeneFromGTF.TranscriptFromGTF tx = gene.addTranscript(transcriptName, transcriptionStart, transcriptionEnd, codingStart, codingEnd, exons.size(), transcriptName, transcriptID, transcriptType);

        for (int i = 0; i < exons.size(); ++i) {
            Exon e = exons.get(i);
            if (e.start > e.end)
				throw new AnnotationException("Exon has 0 or negative extent for " + transcriptDescription);
            if (i > 0 && exons.get(i-1).end >= exons.get(i).start)
				throw new AnnotationException("Exons overlap for " + transcriptDescription);
            tx.addExon(e.start, e.end);
        }


        return tx;


    }


    private void validateGTFRecord(final GTFRecord gtfRecord, final GeneFromGTF gene) {
            if (gene.isPositiveStrand()==gtfRecord.isNegativeStrand())
				throw new AnnotationException("Strand disagreement in GTF file for gene " + gene.getName());
            if (!gene.getContig().equals(gtfRecord.getChromosome()))
				throw new AnnotationException("Chromosome disagreement(" + gene.getContig() + " != " + gtfRecord.getChromosome() +
                        ") in GTF file for gene " + gene.getName());
            if (gtfRecord.getFeatureType().equals(GTFParser.GTFFeature.gene.name()))
				if (gtfRecord.getStart() != gene.getStart() ||
                        gtfRecord.getEnd() != gene.getEnd())
					throw new AnnotationException(String.format("gene GTFRecord(%s) != GeneFromGTF(%s)", gtfRecord.toString(), gene.toString()));

        }

    private Map<String, Collection<GTFRecord>> gatherByGeneName(final Iterator<GTFRecord> gtfRecords) {
        return CollectionUtil.partition(CollectionUtil.makeCollection(gtfRecords),
                new CollectionUtil.Partitioner<GTFRecord, String>() {
                    @Override
                    public String getPartition(final GTFRecord gtfRecord) {
                        return gtfRecord.getGeneName();
                    }
                });
    }

    private Map<Integer, Collection<GTFRecord>> gatherByGeneVersion(final Collection<GTFRecord> gtfRecords) {
        return CollectionUtil.partition(gtfRecords,
                new CollectionUtil.Partitioner<GTFRecord, Integer>() {
                    @Override
                    public Integer getPartition(final GTFRecord gtfRecord) {
                        if (gtfRecord.getGeneVersion() != null)
							return gtfRecord.getGeneVersion();
						else
							return Integer.MIN_VALUE;
                    }
                }
        );
    }

    private Map<String, Collection<GTFRecord>> gatherByTranscriptId(final Collection<GTFRecord> gtfRecords) {
        return CollectionUtil.partition(gtfRecords,
                new CollectionUtil.Partitioner<GTFRecord, String>() {
                    @Override
                    public String getPartition(final GTFRecord gtfRecord) {
                        if (gtfRecord.getTranscriptID() == null)
							throw new RuntimeException("GTFRecord does not have transcriptID: " + gtfRecord);
                        return gtfRecord.getTranscriptID();
                    }
                });
    }

    private static class Exon implements Comparable<Exon> {
        public final int start;
        public final int end;

        public Exon(final int start, final int end) {
            this.start = start;
            this.end = end;
        }

        @Override
        public int compareTo(final Exon o) {
            int ret = Integer.compare(this.start, o.start);
            if (ret != 0)
				return ret;
            return Integer.compare(this.end, o.end);
        }
    }

    private static class GeneAnnotationFilter extends FilteredIterator<GTFRecord> {
        private GeneAnnotationFilter(final Iterator<GTFRecord> underlyingIterator) {
            super(underlyingIterator);
        }

        @Override
        public boolean filterOut(final GTFRecord rec) {
            return GTFParser.GTFFeature.gene.name().equals(rec.getFeatureType());
        }
    }
}
