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
package org.broadinstitute.dropseqrna.cluster;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.RuntimeIOException;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeaderCodec;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeaderMerger;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.yaml.snakeyaml.Yaml;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Argument;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.regex.Pattern;

@CommandLineProgramProperties(
        summary = "To be invoked by the clustering workflow.  " +
                "Read a YAML containing a list of DGEs and some configuration parameters, \n" +
                "and produce a merged, filtered Matrix Market sparse DGE, with genes and cell barcodes in the MM header. \n" +
                "Genes are alphabetically sorted.  Cell barcodes are ordered according to the location of the input DGE \n" +
                "in the YAML, and within a single input DGE are ordered by size (descending).",
        oneLineSummary = "Merge DGEs into a sparse Matrix Market DGE",
        programGroup = DropSeq.class
)
public class MergeDgeSparse
        extends CommandLineProgram {

    @Argument(shortName = "Y", doc="Input file containing list of data sets as for clustering workflow.\n" +
            "The file is expected to contain a 'datasets' list.  Each element of the list will contain:\n\n" +
            "path: the location of the DGE. (required)\n" +
            "name: a prefix to prepend to each cell barcode (plus underscore).  (optional)\n" +
            "cell_count: If present and non-zero, take this many cells from the input DGE, selecting the largest. (optional)\n\n" +
            "YAML argument is")
    public File YAML;

    @Argument(optional = true,
    doc="If set, write a table of cell barcodes and number of transcripts to a possibly gzipped output file.")
    public File CELL_SIZE_OUTPUT_FILE;

    @Argument(optional = true,
    doc = "If set, write a Matrix Market sparse matrix containing the merge of all the DGEs listed in the YAML.  " +
            "At least one of RAW_DGE_OUTPUT_FILE and SCALED_DGE_OUTPUT_FILE must be set.")
    public File RAW_DGE_OUTPUT_FILE;

    @Argument(optional = true,
            doc = "If set, write a Matrix Market sparse matrix containing the merge of all the DGEs listed in the YAML.  " +
                    "The expression for each {gene,cell} is divided by the total expression of the cell.  " +
                    "At least one of RAW_DGE_OUTPUT_FILE and SCALED_DGE_OUTPUT_FILE must be set.")
    public File SCALED_DGE_OUTPUT_FILE;

    @Argument(optional = true,
            doc="If set, write a DGE header that is the result of merging the headers of the input DGEs.")
    public File DGE_HEADER_OUTPUT_FILE;

    @Argument(doc="Remove genes with fewer than this many cells.  This filtering is done after all the cell-based filters.")
    public int MIN_CELLS = Defaults.MIN_CELLS;

    @Argument(doc="Remove cells with fewer than this many genes.")
    public int MIN_GENES = Defaults.MIN_GENES;

    @Argument(doc="Remove cells with fewer than this many transcripts.")
    public int MIN_TRANSCRIPTS = Defaults.MIN_TRANSCRIPTS;

    @Argument(doc="Genes that match one of these regular expressions will be removed.")
    public List<String> FILTERED_GENE_RE;

    @Argument(doc="If specified, files containing lists of cell barcodes, one per line.  " +
            "May be gzipped.  Lines starting with # are ignored.  " +
            "Only the cells in listed in these file(s) are included in the output.  If no files are specified, all cell" +
            "barcodes are included, subject to the other filters.")
    public List<File> SELECTED_CELL_BARCODES;


    @Argument(doc="Controls stringency of DGE header merging.  Only relevant if DGE_HEADER_OUTPUT_FILE is set.")
    public DgeHeaderMerger.Stringency HEADER_STRINGENCY = DgeHeaderMerger.Stringency.STRICT;

    private static final Log LOG = Log.getInstance(MergeDgeSparse.class);

    // Yaml keys organized hierarchically
    static class YamlKeys {
        static final String DATASETS_KEY = "datasets";
        static class DatasetsKeys {
            static final String PATH_KEY = "path";
            static final String CELL_COUNT_KEY = "cell_count";
            static final String NAME_KEY = "name";
        }
    }

    static class Defaults {
        static final int MIN_CELLS = 0;
        static final int MIN_GENES = 400;
        static final int MIN_TRANSCRIPTS = 0;
    }

    public static void main(final String[] args) {
        new MergeDgeSparse().instanceMainWithExit(args);
    }

    private GeneEnumerator geneEnumerator;

    private Set<String> selectedCells = null;

    @Override
    protected int doWork() {
        final Yaml yaml = new Yaml();
        final Map yamlMap = (Map)yaml.load(IOUtil.openFileForReading(YAML));
        geneEnumerator = new GeneEnumerator(FILTERED_GENE_RE);

        @SuppressWarnings("unchecked")
        final List<Map> dataSets = (List<Map>) getRequiredValue(yamlMap, YamlKeys.DATASETS_KEY);

        if (SELECTED_CELL_BARCODES != null && SELECTED_CELL_BARCODES.size() > 0) {
            selectedCells = loadSelectedCellsLists(SELECTED_CELL_BARCODES);
        }

        List<SparseDge> dges = loadDataSets(dataSets);

        writeCellSizesFile(dges);

        writeDgeHeader(dataSets);

        final GeneFiltererSorter geneFiltererSorter = new GeneFiltererSorter(MIN_CELLS, dges);

        int totalNumElements = 0;
        final List<String> cellBarcodes = new ArrayList<>();
        for (final SparseDge dge : dges) {
            totalNumElements += dge.getNumNonZeroEntries();
            for (int i = 0; i < dge.getNumCells(); ++i) {
                cellBarcodes.add(dge.getCellBarcode(i));
            }
        }

        final MergeDgeOutputWriter writer = new MergeDgeOutputWriter(RAW_DGE_OUTPUT_FILE, SCALED_DGE_OUTPUT_FILE,
                totalNumElements, geneFiltererSorter.getSortedGeneNames(), cellBarcodes);

        int cellIndexOffset = 0;
        int numFilteredElements = 0;
        for (final SparseDge dge : dges) {
            LOG.info("Merging DGE " + dge.getFile().getAbsolutePath());
            for (final SparseDge.Triplet triplet: dge.getTriplets()) {
                final int geneIndex = geneFiltererSorter.getOutputGeneIndex(triplet.geneIndex);
                if (geneIndex < 0) {
                    // Gene was filtered by a threshold.
                    ++numFilteredElements;
                    continue;
                }
                final int cellIndex = cellIndexOffset + triplet.cellIndex;
                final double scaled = triplet.value/(double)dge.getNumTranscripts(triplet.cellIndex);
                writer.writeValue(geneIndex, cellIndex, triplet.value, scaled);
            }
            cellIndexOffset += dge.getNumCells();
        }
        writer.close();
        LOG.info(numFilteredElements + " filtered by a gene filter.");
        return 0;
    }

    private Set<String> loadSelectedCellsLists(final List<File> files) {
        final Set<String> ret = new HashSet<>();
        final Pattern comment = Pattern.compile("#");
        final Pattern whitespace = Pattern.compile("\\s");
        for (final File file : files) {
            final BufferedReader reader = IOUtil.openFileForBufferedReading(file);
            String line;
            try {
                while ((line = reader.readLine()) != null) {
                    // Remove trailing comments
                    String[] fields = comment.split(line, 2);
                    if (!fields[0].isEmpty()) {
                        // Remove trailing whitespace
                        fields = whitespace.split(fields[0], 2);
                        if (!fields[0].isEmpty()) {
                            ret.add(fields[0]);
                        }
                    }
                }
            } catch (IOException e) {
                throw new RuntimeIOException("Exception reading " + file.getAbsolutePath(), e);
            }
        }
        return ret;
    }

    @Override
    protected String[] customCommandLineValidation() {
        String[] superMessages = super.customCommandLineValidation();
        if (RAW_DGE_OUTPUT_FILE == null && SCALED_DGE_OUTPUT_FILE == null) {
            return CustomCommandLineValidationHelper.makeValue(superMessages,
                    Collections.singletonList("At least one of RAW_DGE_OUTPUT_FILE and SCALED_DGE_OUTPUT_FILE should be set"));
        } else {
            return superMessages;
        }
    }

    private class GeneFiltererSorter {
        private int numOutputGenes = 0;
        private final int[] geneIdMapping;
        private final List<String> sortedGeneNames;

        public GeneFiltererSorter(final int minCellsPerGene, final List<SparseDge> dges) {
            final int[] cellsPerGene = countCellsPerGene(dges);
            geneIdMapping = new int[geneEnumerator.getGenes().size()];
            final Map<String, Integer> geneMap = new TreeMap<>();
            for (int i = 0; i < geneIdMapping.length; ++i) {
                if (cellsPerGene[i] >= minCellsPerGene) {
                    geneMap.put(geneEnumerator.getGeneName(i), i);
                } else {
                    geneIdMapping[i] = -1;
                }
            }
            for (final Map.Entry<String, Integer> entry : geneMap.entrySet()) {
                geneIdMapping[entry.getValue()] = numOutputGenes++;
            }
            sortedGeneNames = Collections.unmodifiableList(new ArrayList<String>(geneMap.keySet()));
        }

        private int[] countCellsPerGene(final List<SparseDge> dges) {
            final int[] cellsPerGene = new int[geneEnumerator.getGenes().size()];
            for (final SparseDge dge : dges) {
                for (final SparseDge.Triplet triplet : dge.getTriplets()) {
                    ++cellsPerGene[triplet.geneIndex];
                }
            }
            return cellsPerGene;
        }

        public int getNumOutputGenes() {
            return numOutputGenes;
        }

        public int getOutputGeneIndex(final int originalGeneIndex) {
            return geneIdMapping[originalGeneIndex];
        }

        public List<String> getSortedGeneNames() {
            return sortedGeneNames;
        }
    }

    private void writeDgeHeader(final List<Map> dataSets) {
        if (DGE_HEADER_OUTPUT_FILE != null) {
            LOG.info("Writing " + DGE_HEADER_OUTPUT_FILE.getAbsolutePath());
            final List<File> inputDges = new ArrayList<>(dataSets.size());
            final List<String> prefixes = new ArrayList<>(dataSets.size());
            for (final Map dataSet : dataSets) {
                inputDges.add(new File((String)dataSet.get(YamlKeys.DatasetsKeys.PATH_KEY)));
                prefixes.add((String)getValueOrDefault(dataSet, YamlKeys.DatasetsKeys.NAME_KEY, ""));
            }
            new DgeHeaderCodec().encode(DGE_HEADER_OUTPUT_FILE, DgeHeaderMerger.mergeDgeHeaders(inputDges, prefixes, HEADER_STRINGENCY));
        }
    }

    private void writeCellSizesFile(final List<SparseDge> dges) {
        if (CELL_SIZE_OUTPUT_FILE != null) {
            LOG.info("Writing " + CELL_SIZE_OUTPUT_FILE.getAbsolutePath());
            final CellSizeWriter cellSizeWriter = new CellSizeWriter(CELL_SIZE_OUTPUT_FILE);
            for (final SparseDge dge : dges) {
                for (int i = 0; i < dge.getNumCells(); ++i) {
                    cellSizeWriter.writeSize(dge.getCellBarcode(i), dge.getNumTranscripts(i));
                }
            }
            cellSizeWriter.close();
        }
    }

    private List<SparseDge> loadDataSets(final List datasets) {
        if (datasets.size() == 0) {
            throw new RuntimeException("List of datasets is empty in " + YAML.getAbsolutePath());
        }
        final ArrayList<SparseDge> dges = new ArrayList<>(datasets.size());

        for (final Object dataset : datasets) {
            if (! (dataset instanceof Map)) {
                throw new RuntimeException("dataset element is not map in " + YAML.getAbsolutePath());
            }
            final Map datasetMap = (Map)dataset;
            final File dgePath = new File((String)getRequiredValue(datasetMap, YamlKeys.DatasetsKeys.PATH_KEY));
            String prefix = (String)getValueOrDefault(datasetMap, YamlKeys.DatasetsKeys.NAME_KEY, "");
            final SparseDge dge = new SparseDge(dgePath, geneEnumerator);
            LOG.info(String.format("Loaded %d cells from %s", dge.getNumCells(), dgePath.getAbsolutePath()));

            if (!prefix.isEmpty()) {
                dge.prefixCellBarcodes(prefix + "_");
            }
            Integer cell_count = (Integer)getValueOrDefault(datasetMap, YamlKeys.DatasetsKeys.CELL_COUNT_KEY, 0);
            if (cell_count > 0 && cell_count < dge.getNumCells()) {
                dge.discardSmallestCells(cell_count);
            }

            if (selectedCells != null) {
                dge.retainOnlyTheseCells(selectedCells);
            }

            if (MIN_GENES > 0) {
                int numCells = dge.getNumCells();
                dge.discardCellsWithFewGenes(MIN_GENES);
                LOG.info(String.format("Discarded %d cells with fewer than %d genes from %s", numCells - dge.getNumCells(), MIN_GENES, dgePath.getAbsolutePath()));
            }
            if (MIN_TRANSCRIPTS > 0) {
                int numCells = dge.getNumCells();
                dge.discardCellsWithFewTranscripts(MIN_TRANSCRIPTS);
                LOG.info(String.format("Discarded %d cells with fewer than %d transcripts from %s", numCells - dge.getNumCells(), MIN_TRANSCRIPTS, dgePath.getAbsolutePath()));
            }
            dges.add(dge);
        }
        return dges;
    }

    private Object getRequiredValue(final Map map, final String key) {
        final Object ret = map.get(key);
        if (ret == null) {
            throw new RuntimeException(YAML.getAbsolutePath() + " does not contain key " + key);
        }
        return ret;
    }

    private Object getValueOrDefault(final Map map, final String key, final Object defaultValue) {
        final Object ret = map.get(key);
        if (ret == null) {
            return defaultValue;
        } else {
            return ret;
        }
    }
}
