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

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.*;
import picard.annotation.AnnotationException;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.File;
import java.util.List;
import java.util.Map;

/**
 * Simple class that merely parses a GTF file into GTFRecords.
 * All the smarts about gathering lines into genes, filtering, etc., is done elsewhere
 */
public class GTFParser extends IterableOnceIterator<GTFRecord> {

    // These are in the order that columns appear in GTF format.
    public enum GTFFlatColumns{CHROMOSOME, SOURCE, FEATURE, START, END, SCORE, STRAND, FRAME, ATTRIBUTE}

    /**
     * Not a comprehensive list, just the common ones
     */
    public enum GTFFeature{CDS, UTR, exon, gene, start_codon, stop_codon, transcript}

    private static final Log LOG = Log.getInstance(GTFParser.class);

    private static final String[] GTFColumnLabels = new String[GTFFlatColumns.values().length];

    static {
        for (int i = 0; i < GTFColumnLabels.length; ++i) {
            GTFColumnLabels[i] = GTFFlatColumns.values()[i].name();
        }
    }

    private final File gtfFile;
    private final CloseableIterator<TabbedTextFileWithHeaderParser.Row> it;
    private final ProgressLogger progressLogger;
    private final ValidationStringency validationStringency;

    public GTFParser(final File gtfFile, final ValidationStringency validationStringency) {
        this.gtfFile = gtfFile;
        this.it = new TabbedTextFileWithHeaderParser(gtfFile, GTFColumnLabels).iterator();
        progressLogger = new ProgressLogger(LOG, 100000, "read", "GTF records");
        this.validationStringency = validationStringency;
    }

    @Override
    public void close() {
        it.close();
    }

    @Override
    public boolean hasNext() {
        return it.hasNext();
    }

    @Override
    public GTFRecord next() {
        final TabbedTextFileWithHeaderParser.Row row = it.next();
        if (row.getFields().length != GTFColumnLabels.length) {
            throw new AnnotationException("Wrong number of fields in GTF file " + gtfFile + " at line " +
                    row.getCurrentLine());
        }
        final GTFRecord ret = parseLine(row);
        if (validationStringency != ValidationStringency.SILENT) {
            final List<String> errors = ret.validate();
            if (errors != null && !errors.isEmpty()) {
                final String message = String.format(
                        "Invalid GTF line: \n%s\nProblems:\n%s",
                        row.getCurrentLine(),
                        CollectionUtil.join(errors, "\n"));
                if (validationStringency == ValidationStringency.STRICT) {
                    throw new AnnotationException(message);
                } else {
                    LOG.warn(message);
                }
            }
        }
        progressLogger.record(ret.getChromosome(), ret.getStart());
        return ret;
    }

    @Override
    public void remove() {
        it.remove();
    }

    private GTFRecord parseLine (TabbedTextFileWithHeaderParser.Row row) {
        final String attributes= row.getField(GTFFlatColumns.ATTRIBUTE.name());
        final Map<String, String> attributesMap = AnnotationUtils.getInstance().parseOptionalFields(attributes);

        final String chromosome = row.getField(GTFFlatColumns.CHROMOSOME.name());
        final int start = Integer.parseInt(row.getField(GTFFlatColumns.START.name()));
        final int end = Integer.parseInt(row.getField(GTFFlatColumns.END.name()));
        final String strand = row.getField(GTFFlatColumns.STRAND.name());
        final String featureType = row.getField(GTFFlatColumns.FEATURE.name());

        final String geneName=attributesMap.get("gene_name");
        final String geneID=attributesMap.get("gene_id");
        final String transcriptName=attributesMap.get("transcript_name");
        final String transcriptID=attributesMap.get("transcript_id");
        final String transcriptType=attributesMap.get("gene_biotype");
        final String geneVersionString = attributesMap.get("gene_version");
        final Integer geneVersion;
        if (geneVersionString != null) {
            geneVersion = new Integer(geneVersionString);
        } else {
            geneVersion = null;
        }

        final boolean negativeStrand=strand.equals("-");

        return new GTFRecord (chromosome, start, end, negativeStrand, geneID, geneName, transcriptName, transcriptID,
                transcriptType, featureType, geneVersion);
    }
}
