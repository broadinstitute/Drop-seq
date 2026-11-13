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
package org.broadinstitute.dropseqrna.barnyard;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.metrics.BamTagHistogram;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.readiterators.GeneFunctionProcessor;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.RuntimeIOException;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        summary = "Determine the cell barcodes that have a minimum number of transcripts.  For a multi-species BAM, " +
            "the minimum must be reached in transcripts for a single species.",
        oneLineSummary = "Determine the cell barcodes that have a minimum number of transcripts.",
        programGroup = DropSeq.class)
public class SelectCellsByNumTranscripts
        extends GeneFunctionCommandLineBase {

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze. This argument can accept wildcards, or a file with the suffix .bam_list that contains the locations of multiple BAM files", minElements = 1)
    public List<File> INPUT;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="List of cell barcodes, one per line")
    public File OUTPUT;

    @Argument(shortName = StandardOptionDefinitions.METRICS_FILE_SHORT_NAME, optional = true,
            doc="If specified, produce summary metrics")
    public File METRICS;

    @Argument(doc="Select cells with at least this many transcripts")
    public Integer MIN_TRANSCRIPTS_PER_CELL;

    @Argument(doc="Limit cells to those with at least this many reads.", optional = true)
    public Integer MIN_READS_PER_CELL;

    @Argument(doc="If specified, transcripts are counted per species, and the MIN_TRANSCRIPTS_PER_CELL threshold must be reached by transcripts for a single species before a cell is selected.",
            optional = true)
    public List<String> ORGANISM;

    @Argument(doc="If set, cells with minimum number of reads are written to this file.", optional = true)
    public File OUTPUT_INTERIM_CELLS;

    @Argument(doc="If set, read cells from this file rather than filtering BAM for cells with minimum number of reads.", optional = true)
    public File INPUT_INTERIM_CELLS;

    @Argument(doc="The cell barcode tag.  If there are no reads with this tag, the program will assume that all reads belong to the same cell and process in single sample mode.")
    public String CELL_BARCODE_TAG="XC";

    @Argument(doc="The molecular barcode tag.")
    public String MOLECULAR_BARCODE_TAG="XM";

    @Argument(doc="The map quality of the read to be included.")
    public int READ_MQ=10;

    @Argument(doc="The edit distance that molecular barcodes should be combined at within a gene.")
    public Integer EDIT_DISTANCE=1;

    @Argument(doc = "Number of threads used to read and decode the input BAMs, for both the read count pass and the UMI pass.  " +
            "More threads help most with several input BAMs, roughly one thread per BAM.  A single BAM speeds up as threads are added, up to about 2 threads.")
    public int NUM_THREADS = 1;

    private static final Log log = Log.getInstance(SelectCellsByNumTranscripts.class);
    static final String ORGANISM_SEPARATOR = "::";

    @Override
    protected int doWork() {
    	this.INPUT = FileListParsingUtils.expandFileList(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        if (MIN_READS_PER_CELL == null || MIN_READS_PER_CELL < MIN_TRANSCRIPTS_PER_CELL)
			MIN_READS_PER_CELL = MIN_TRANSCRIPTS_PER_CELL;
        // Cells with enough reads are found while the reads are being sorted, unless the list is given.
        List<String> cellBarcodes = INPUT_INTERIM_CELLS == null ? null : readBarcodes(INPUT_INTERIM_CELLS);

        SamHeaderAndIterator headerAndIterator = wrapWithReadCounter(SamFileMergeUtil.mergeInputPaths(
                IOUtil.filesToPaths(this.INPUT), false, SamReaderFactory.makeDefault(), this.NUM_THREADS));

        if (ORGANISM != null && !ORGANISM.isEmpty())
            headerAndIterator = new SamHeaderAndIterator(headerAndIterator.header, new PrefixGeneWithOrganismIterator(headerAndIterator.iterator));

        // gene/exon tags are sorted first, followed by cells
        UMIIterator umiIterator = new UMIIterator.UMIIteratorBuilder(headerAndIterator, GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG,
        		this.STRAND_STRATEGY, this.LOCUS_FUNCTION_LIST, this.FUNCTIONAL_STRATEGY, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG,
        		this.READ_MQ).setCellBarcodes(cellBarcodes).build();

        cellBarcodes = selectCellsByReadCount(cellBarcodes);

        // Per organism, the transcripts and genes of each cell.  Groups of cells that are not in the list are not counted.
        final String[] genePrefixes = getGenePrefixes();
        final List<Map<String, DigitalExpression.DESummary>> summaries = new ArrayList<>();
        for (int i = 0; i < Math.max(1, genePrefixes.length); i++)
            summaries.add(DigitalExpression.initializeSummary(cellBarcodes));

        UMICollection batch;
        while ((batch=umiIterator.next())!=null) {
            if (batch.isEmpty())
				continue;
            final int organism = getOrganismIndex(genePrefixes, batch.getGeneName());
            final DigitalExpression.DESummary summary = organism < 0 ? null : summaries.get(organism).get(batch.getCellBarcode());
            if (summary == null)
				continue;
            // every group is one gene in one cell, so it can be counted as it comes
            summary.NUM_TRANSCRIPTS += batch.getDigitalExpression(0, this.EDIT_DISTANCE, false);
            summary.NUM_GENES++;
        }

        final Map<String, Integer> transcriptsPerCell = getTranscriptsPerSelectedCell(summaries);

        log.info("Found " + transcriptsPerCell.size() + " cells with enough transcripts");

        final List<Map.Entry<String, Integer>> sortedCells = new ArrayList<>(transcriptsPerCell.entrySet());
        sortedCells.sort(Map.Entry.comparingByValue(Comparator.reverseOrder()));

        final List<String> finalBarcodes = new ArrayList<>(sortedCells.size());
        for (final Map.Entry<String, Integer> entry: sortedCells)
			finalBarcodes.add(entry.getKey());

        writeBarcodes(OUTPUT, finalBarcodes);
        log.info("Wrote cell barcodes to " + OUTPUT.getAbsolutePath());
        CloserUtil.close(umiIterator);

        if (METRICS != null) {
            final Metrics m = new Metrics();
            for (final Map<String, DigitalExpression.DESummary> organismSummaries: summaries)
				for (final DigitalExpression.DESummary summary: organismSummaries.values())
					m.accumulate(summary, transcriptsPerCell.containsKey(summary.CELL_BARCODE));
            MetricsFile<Metrics, Integer> out = getMetricsFile();
            out.addMetric(m);
            out.write(METRICS);
        }
        return 0;
    }


    /** Counts the reads of each cell as they pass, in the same way as {@link org.broadinstitute.dropseqrna.metrics.BamTagHistogram}. */
    private class ReadCountingIterator
    implements CloseableIterator<SAMRecord> {
        private final CloseableIterator<SAMRecord> it;
        private final ObjectCounter<String> counts = new ObjectCounter<>();
        private boolean exhausted = false;

        public ReadCountingIterator(final CloseableIterator<SAMRecord> it) {
            this.it = it;
        }

        public ObjectCounter<String> getCounts() { return counts; }

        /** @return true once all reads have been passed on.  The underlying iterator may be closed by then. */
        public boolean isExhausted() { return exhausted; }

        @Override
        public void close() { it.close(); }

        @Override
        public boolean hasNext() {
            final boolean hasNext = it.hasNext();
            if (!hasNext)
                exhausted = true;
            return hasNext;
        }

        @Override
        public SAMRecord next() {
            final SAMRecord rec = it.next();
            BamTagHistogram.countRecord(counts, rec, CELL_BARCODE_TAG, READ_MQ, false);
            return rec;
        }

        @Override
        public void remove() { it.remove(); }
    }

    private class PrefixGeneWithOrganismIterator
    implements CloseableIterator<SAMRecord> {
        private final CloseableIterator<SAMRecord> it;

        private final String[] referenceSequencePrefix = new String[ORGANISM.size()];
        private final String[] genePrefix = new String[ORGANISM.size()];

        public PrefixGeneWithOrganismIterator(final CloseableIterator<SAMRecord> it) {
            this.it = it;
            for (int i = 0; i < ORGANISM.size(); ++i) {
                referenceSequencePrefix[i] = ORGANISM.get(i) + "_";
                genePrefix[i] = ORGANISM.get(i) + ORGANISM_SEPARATOR;
            }
        }

        @Override
        public void close() { it.close(); }

        @Override
        public boolean hasNext() { return it.hasNext(); }

        @Override
        public SAMRecord next() {
            final SAMRecord rec = it.next();
            final String reference = rec.getReferenceName();
            for (int i = 0; i < referenceSequencePrefix.length; ++i)
				if (reference.startsWith(referenceSequencePrefix[i])) {
                    final String geneExon = rec.getStringAttribute(GENE_NAME_TAG);
                    if (geneExon != null) {
                        // A read that overlaps several genes has a comma separated list, every one needs the prefix.
                        rec.setAttribute(GENE_NAME_TAG, genePrefix[i] +
                                geneExon.replace(GeneFunctionProcessor.DELIMITER, GeneFunctionProcessor.DELIMITER + genePrefix[i]));
                        break;
                    }
                }
            return rec;
        }

        @Override
        public void remove() { it.remove(); }
    }

    public static void writeBarcodes(final File file, final List<String> barcodes) {
        final BufferedWriter writer = IOUtil.openFileForBufferedWriting(file);
        try {
            for (final String barcode : barcodes) {
                writer.write(barcode);
                writer.newLine();
            }
            writer.close();
        } catch (IOException e) {
            throw new RuntimeIOException("Exception writing " + file.getAbsolutePath(), e);
        }
    }

    public static List<String> readBarcodes(final File file) {
        try {
            IOUtil.assertFileIsReadable(file);
            final BufferedReader reader = IOUtil.openFileForBufferedReading(file);
            final List<String> ret = new ArrayList<>();
            String barcode;
            while ((barcode = reader.readLine()) != null)
				if (!barcode.isEmpty())
					ret.add(barcode);
            CloserUtil.close(reader);
            return ret;
        } catch (IOException e) {
            throw new RuntimeIOException("Exception reading " + file.getAbsolutePath());
        }
    }

    /** The prefix that the genes of each organism are given, empty if there is a single organism. */
    private String[] getGenePrefixes() {
        if (ORGANISM == null)
			return new String[0];
        return ORGANISM.stream().map(organism -> organism + ORGANISM_SEPARATOR).toArray(String[]::new);
    }

    /** @return the index of the organism of the gene, or -1 if the gene is not an organism's.  With a single organism, it is always 0. */
    private static int getOrganismIndex(final String[] genePrefixes, final String gene) {
        if (genePrefixes.length == 0)
			return 0;
        for (int i = 0; i < genePrefixes.length; ++i)
			if (gene.startsWith(genePrefixes[i]))
				return i;
        return -1;
    }

    /**
     * The cells with at least MIN_TRANSCRIPTS_PER_CELL transcripts in an organism, with their most transcripts in any organism.
     * The cells are in the order of the iteration of the summaries, which decides the order of cells with the same number of transcripts.
     */
    private Map<String, Integer> getTranscriptsPerSelectedCell(final List<Map<String, DigitalExpression.DESummary>> summaries) {
        final Map<String, Integer> result = new HashMap<>();
        for (final DigitalExpression.DESummary first: summaries.get(0).values()) {
            int transcripts = 0;
            for (final Map<String, DigitalExpression.DESummary> organismSummaries: summaries)
				transcripts = Math.max(transcripts, organismSummaries.get(first.CELL_BARCODE).NUM_TRANSCRIPTS);
            if (transcripts >= MIN_TRANSCRIPTS_PER_CELL)
				result.put(first.CELL_BARCODE, transcripts);
        }
        return result;
    }

    /** Counts the reads of each cell as they pass, null if the list of cells was given. */
    private ReadCountingIterator readCounter = null;

    /**
     * If the cells with enough reads are not given, count the reads of each cell as they pass through, so that
     * the input does not have to be read an extra time to find the cells.
     */
    private SamHeaderAndIterator wrapWithReadCounter(final SamHeaderAndIterator headerAndIterator) {
        if (INPUT_INTERIM_CELLS != null)
            return headerAndIterator;
        readCounter = new ReadCountingIterator(headerAndIterator.iterator);
        return new SamHeaderAndIterator(headerAndIterator.header, readCounter);
    }

    /**
     * The cells with at least MIN_READS_PER_CELL reads, and write them to OUTPUT_INTERIM_CELLS if that is set.
     * Call after all of the reads have passed through the read counter.
     * @param givenCellBarcodes the cells from INPUT_INTERIM_CELLS, null if the reads were counted.
     */
    private List<String> selectCellsByReadCount(final List<String> givenCellBarcodes) {
        final List<String> cellBarcodes;
        if (readCounter == null)
            cellBarcodes = givenCellBarcodes;
        else {
            if (!readCounter.isExhausted())
                throw new IllegalStateException("The reads were not all read before the cells were selected.");
            cellBarcodes = new BarcodeListRetrieval().getCoreBarcodesByReadCount(readCounter.getCounts(), MIN_READS_PER_CELL);
        }
        if (OUTPUT_INTERIM_CELLS != null)
            writeBarcodes(OUTPUT_INTERIM_CELLS, cellBarcodes);
        return cellBarcodes;
    }

    @Override
    protected String[] customCommandLineValidation() {
        List<String> thisErrors = null;
        if (NUM_THREADS < 1)
			thisErrors = Arrays.asList("NUM_THREADS must be >= 1");
        else if (ORGANISM.size() != new HashSet<>(ORGANISM).size())
			thisErrors = Arrays.asList("Duplicates not allow in ORGANISM argument");
		else
			for (final String organism : ORGANISM)
				if (organism.contains(ORGANISM_SEPARATOR))
					thisErrors = Arrays.asList("'" + ORGANISM_SEPARATOR + "' not allowed in ORGANISM argument");
        return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), thisErrors);
    }

    public static void main(final String[] args) {
        new SelectCellsByNumTranscripts().instanceMainWithExit(args);
    }

    public static class Metrics
    extends MetricBase {
        public int NUM_TRANSCRIPTS;
        public int NUM_GENES;
        public int NUM_TRANSCRIPTS_SELECTED_CELLS;
        public int NUM_GENES_SELECTED_CELLS;

        public void accumulate(final DigitalExpression.DESummary summary, boolean isSelected) {
            NUM_TRANSCRIPTS += summary.NUM_TRANSCRIPTS;
            NUM_GENES += summary.NUM_GENES;
            if (isSelected) {
                NUM_TRANSCRIPTS_SELECTED_CELLS += summary.NUM_TRANSCRIPTS;
                NUM_GENES_SELECTED_CELLS += summary.NUM_GENES;
            }
        }

        public void accumulate(final Metrics other) {
            NUM_TRANSCRIPTS += other.NUM_TRANSCRIPTS;
            NUM_GENES += other.NUM_GENES;
            NUM_TRANSCRIPTS_SELECTED_CELLS += other.NUM_TRANSCRIPTS_SELECTED_CELLS;
            NUM_GENES_SELECTED_CELLS += other.NUM_GENES_SELECTED_CELLS;
        }
    }
}
