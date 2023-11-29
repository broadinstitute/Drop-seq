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
import java.util.Set;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator;

import htsjdk.samtools.SAMRecord;
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

    private static final Log log = Log.getInstance(SelectCellsByNumTranscripts.class);
    static final String ORGANISM_SEPARATOR = "::";

    @Override
    protected int doWork() {
    	this.INPUT = FileListParsingUtils.expandFileList(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        if (MIN_READS_PER_CELL == null || MIN_READS_PER_CELL < MIN_TRANSCRIPTS_PER_CELL)
			MIN_READS_PER_CELL = MIN_TRANSCRIPTS_PER_CELL;
        final List<String> cellBarcodes;
        if (INPUT_INTERIM_CELLS == null)
			cellBarcodes = new BarcodeListRetrieval().getListCellBarcodesByReadCount(
                    this.INPUT, this.CELL_BARCODE_TAG, this.READ_MQ, this.MIN_READS_PER_CELL, null);
		else
			cellBarcodes = readBarcodes(INPUT_INTERIM_CELLS);

        if (OUTPUT_INTERIM_CELLS != null)
			writeBarcodes(OUTPUT_INTERIM_CELLS, cellBarcodes);

        SamHeaderAndIterator headerAndIterator = SamFileMergeUtil.mergeInputs(this.INPUT, false);

        final MapContainer mapContainer;
        if (ORGANISM != null && !ORGANISM.isEmpty()) {
            headerAndIterator = new SamHeaderAndIterator(headerAndIterator.header, new PrefixGeneWithOrganismIterator(headerAndIterator.iterator));
            mapContainer = new MultiOrganismMapContainer(cellBarcodes);
        } else
			mapContainer = new SingleOrganismMapContainer(cellBarcodes);

        // gene/exon tags are sorted first, followed by cells
        UMIIterator umiIterator = new UMIIterator(headerAndIterator, GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG,
        		this.STRAND_STRATEGY, this.LOCUS_FUNCTION_LIST, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG,
        		this.READ_MQ, false, cellBarcodes);


        String gene = null;

        UMICollection batch;
        while ((batch=umiIterator.next())!=null) {
            if (batch.isEmpty())
				continue;
            String currentGene = batch.getGeneName();
            // if just starting the loop
            if (gene==null) gene=currentGene;
            // you've gathered all the data for the gene, write it out and start on the next.
            if (!gene.equals(currentGene))
				mapContainer.addToSummary(gene);
            // Accumulate this gene
            gene=currentGene;
            mapContainer.countExpression(gene, batch.getCellBarcode(), batch.getDigitalExpression(0, this.EDIT_DISTANCE, false));

        }
        // write out remainder
        mapContainer.addToSummary(gene);

        final Map<String, Integer> transcriptsPerCell = mapContainer.getTranscriptCountForCellBarcodesOverTranscriptThreshold(MIN_TRANSCRIPTS_PER_CELL);

        log.info("Found " + transcriptsPerCell.size() + " cells with enough transcripts");

        final Map.Entry<String, Integer>[] transcriptsPerCellArray = transcriptsPerCell.entrySet().toArray(new Map.Entry[transcriptsPerCell.size()]);
        Arrays.sort(transcriptsPerCellArray, new EntryComparator());

        final List<String> finalBarcodes = new ArrayList<>(transcriptsPerCellArray.length);
        for (final Map.Entry<String, Integer> entry: transcriptsPerCellArray)
			finalBarcodes.add(entry.getKey());

        writeBarcodes(OUTPUT, finalBarcodes);
        log.info("Wrote cell barcodes to " + OUTPUT.getAbsolutePath());
        CloserUtil.close(umiIterator);

        if (METRICS != null) {
            final Metrics m = mapContainer.accumulateMetrics(transcriptsPerCell.keySet());
            MetricsFile<Metrics, Integer> out = getMetricsFile();
            out.addMetric(m);
            out.write(METRICS);
        }
        return 0;
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
                        rec.setAttribute(GENE_NAME_TAG, genePrefix[i] + geneExon);
                        break;
                    }
                }
            return rec;
        }

        @Override
        public void remove() { it.remove(); }
    }

    private interface MapContainer {
        void addToSummary(final String gene);
        void countExpression(final String gene, final String cellBarcode, final int molBCCount);
        Map<String, Integer> getTranscriptCountForCellBarcodesOverTranscriptThreshold(final int minNumTranscripts);
        Metrics accumulateMetrics(final Set<String> selectedCellBarcodes);
    }

    private class SingleOrganismMapContainer
            implements MapContainer {
        final Map<String, Integer> countMap = new HashMap<>();

        final Map<String, DigitalExpression.DESummary> summaryMap;

        public SingleOrganismMapContainer(final List<String> cellBarcodes) {
            summaryMap = DigitalExpression.initializeSummary(cellBarcodes);
        }

        @Override
        public void addToSummary(final String gene) {
            DigitalExpression.addToSummary(getZeroValueMap(), countMap, summaryMap);
            countMap.clear();
        }

        /**
         * A bit of a work around for DESummary requiring a map of read counts that we don't record here because we don't use it.
         * This creates a map with all cell barcode keys, with all values set to 0.
         * Collections.emptyMap() wasn't cutting it, and changing DESummary to take null values seemed like a bad idea.
         */
        private Map<String, Integer> getZeroValueMap() {
        	Map<String, Integer> result = new HashMap<>();
        	for (String k: countMap.keySet())
				result.put(k, 0);
        	return result;
        }

        @Override
        public void countExpression(final String gene, final String cellBarcode, final int molBCCount) {
            countMap.put(cellBarcode, molBCCount);
        }

        @Override
        public Map<String, Integer> getTranscriptCountForCellBarcodesOverTranscriptThreshold(final int minNumTranscripts) {
            final Map<String, Integer> ret = new HashMap<>();
            for (final DigitalExpression.DESummary summary : summaryMap.values())
				if (summary.NUM_TRANSCRIPTS >= minNumTranscripts)
					ret.put(summary.CELL_BARCODE, summary.NUM_TRANSCRIPTS);
            return ret;
        }

        @Override
        public Metrics accumulateMetrics(Set<String> selectedCellBarcodes) {
            final Metrics m = new Metrics();
            for (final String cellBarcode: summaryMap.keySet()) {
                m.accumulate(summaryMap.get(cellBarcode), selectedCellBarcodes.contains(cellBarcode));
            }
            return m;
        }
    }

    private class MultiOrganismMapContainer
            implements MapContainer {
        final SingleOrganismMapContainer[] innerMapContainer = new SingleOrganismMapContainer[ORGANISM.size()];
        final String[] genePrefixes = new String[ORGANISM.size()];

        public MultiOrganismMapContainer(final List<String> cellBarcodes) {
            for (int i = 0; i < ORGANISM.size(); ++i) {
                innerMapContainer[i] = new SingleOrganismMapContainer(cellBarcodes);
                genePrefixes[i] = ORGANISM.get(i) + ORGANISM_SEPARATOR;
            }
        }

        private int getOrganismIndex(final String gene) {
            for (int i = 0; i < genePrefixes.length; ++i)
				if (gene.startsWith(genePrefixes[i]))
					return i;
            return -1;
        }

        private SingleOrganismMapContainer getInnerMapContainer(final String gene) {
            final int index = getOrganismIndex(gene);
            if (index < 0)
				throw new IllegalArgumentException("Gene '" + gene + "' is not countable");
            return innerMapContainer[index];
        }

        private boolean countableGene(final String gene) {
            return getOrganismIndex(gene) >= 0;
        }

        @Override
        public void addToSummary(final String gene) {
            if (countableGene(gene))
				getInnerMapContainer(gene).addToSummary(gene);
        }

        @Override
        public void countExpression(final String gene, final String cellBarcode, final int molBCCount) {
            if (countableGene(gene))
				getInnerMapContainer(gene).countExpression(gene, cellBarcode, molBCCount);
        }

        @Override
        public Map<String, Integer> getTranscriptCountForCellBarcodesOverTranscriptThreshold(final int minNumTranscripts) {
            Map<String, Integer> ret = null;
            for (SingleOrganismMapContainer somc : innerMapContainer) {
                final Map<String, Integer> organismSet = somc.getTranscriptCountForCellBarcodesOverTranscriptThreshold(minNumTranscripts);
                if (ret == null)
					ret = organismSet;
				else
					for (final Map.Entry<String, Integer> entry : organismSet.entrySet()) {
                        final String cellBarcode = entry.getKey();
                        final Integer numTranscripts = ret.get(cellBarcode);
                        final Integer newNumTranscripts = entry.getValue();
                        if (numTranscripts == null || numTranscripts < newNumTranscripts)
							ret.put(cellBarcode, newNumTranscripts);
                    }
            }
            return ret;
        }

        @Override
        public Metrics accumulateMetrics(Set<String> selectedCellBarcodes) {
            final Metrics m = innerMapContainer[0].accumulateMetrics(selectedCellBarcodes);
            for (int i = 1; i < innerMapContainer.length; i++) {
                m.accumulate(innerMapContainer[i].accumulateMetrics(selectedCellBarcodes));
            }
            return m;
        }
    }

    /**
     * Sorted in descending order by number of transcripts.
     */
    private class EntryComparator implements Comparator<Map.Entry<String, Integer>> {
        @Override
        public int compare(final Map.Entry<String, Integer> o1, final Map.Entry<String, Integer> o2) {
            return o2.getValue().compareTo(o1.getValue());
        }
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

    @Override
    protected String[] customCommandLineValidation() {
        List<String> thisErrors = null;
        if (ORGANISM.size() != new HashSet<>(ORGANISM).size())
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
