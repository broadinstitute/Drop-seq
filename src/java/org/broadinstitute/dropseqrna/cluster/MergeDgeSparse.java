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
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeaderCodec;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeaderMerger;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.matrixmarket.MatrixMarketConstants;
import org.broadinstitute.dropseqrna.matrixmarket.MatrixMarketWriter;
import org.yaml.snakeyaml.Yaml;
import picard.cmdline.CommandLineProgram;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

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
            "cell_count: If present and non-zero, take this many cells from the input DGE, selecting the largest. (optional)" +
            "cell_barcode_path: File containing non-prefixed cell barcodes to include.  (optional)\n\n" +
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

    @Argument(optional = true,
    doc="If set, write a list of cell barcodes that have been filtered by any of the filtering mechanisms.")
    public File DISCARDED_CELLS_FILE;

    @Argument(doc="Remove genes with fewer than this many cells.  This filtering is done after all the cell-based filters." +
            "Must be >= 1.")
    public int MIN_CELLS = Defaults.MIN_CELLS;

    @Argument(doc="Remove cells with fewer than this many genes.  Must be >= 1.")
    public int  MIN_GENES = Defaults.MIN_GENES;

    @Argument(doc="Remove cells with fewer than this many transcripts.  Must be >= 1.")
    public int MIN_TRANSCRIPTS = Defaults.MIN_TRANSCRIPTS;

    @Argument(doc="Genes that match one of these regular expressions will be removed.")
    public List<String> FILTERED_GENE_RE;

    @Argument(doc="If specified, files containing lists of cell barcodes, one per line.  " +
    		"This set of cell barcodes should contain prefixes matching those specified in the YAML file [name attribute]. " +
    		"If you wish to filter with cell barcodes without prefixes, specify these per-experiment in the YAML. " +
            "May be gzipped.  Lines starting with # are ignored.  " +
            "Only the cells in listed in these file(s) are included in the output.  If no files are specified, all cell " +
            "barcodes are included, subject to the other filters.")
    public List<File> CELL_BC_FILE;

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
            static final String CELL_BARCODE_PATH_KEY = "cell_barcode_path";
        }
    }

    static class Defaults {
        static final int MIN_CELLS = 1;
        static final int MIN_GENES = 400;
        static final int MIN_TRANSCRIPTS = 1;
    }

    private static class DgeDescr {
        final File file;
        final String name;
        /** The cell barcodes that have survived thresholds */
        final ArrayList<String> cellBarcodes;
        /** cell barcodes that have been discarded by a threshold. */
        final List<String> discardedCells;
        /** compact representation of gene-cell barcode matrix that have non-zero expression  */
        final DgeBitSet dgeBitSet;
        /** mapping from index of cell barcodes in the input DGE to index in the output (relative to this DGE only) */
        CellBarcodeIndexMap cellBarcodeIndexMap = null;

        public DgeDescr(File file, String name, List<String> cellBarcodes, Collection<String> discardedCells,
                        final DgeBitSet dgeBitSet) {
            this.file = file;
            this.name = name;
            this.cellBarcodes = new ArrayList<>(cellBarcodes);
            this.discardedCells = new ArrayList<>(discardedCells);
            this.dgeBitSet = dgeBitSet;
        }
    }

    public static void main(final String[] args) {
        new MergeDgeSparse().instanceMainWithExit(args);
    }

    private GeneEnumerator geneEnumerator;
    private final CellsPerGeneCounter cellsPerGeneCounter = new CellsPerGeneCounter();


    private Set<String> selectedCells = null;

    @Override
    protected int doWork() {
        final Yaml yaml = new Yaml();
        @SuppressWarnings("unchecked") final Map<String, Object> yamlMap = (Map<String, Object>)yaml.load(IOUtil.openFileForReading(YAML));
        geneEnumerator = new GeneEnumerator(FILTERED_GENE_RE);

        //noinspection unchecked
        final List<Map<String, Object>> dataSets = (List<Map<String, Object>>) getRequiredValue(yamlMap, YamlKeys.DATASETS_KEY);

        if (CELL_BC_FILE != null && CELL_BC_FILE.size() > 0)
            selectedCells = loadSelectedCellsLists(CELL_BC_FILE);

        // Read each DGE, apply MIN_GENES and MIN_TRANSCRIPTS thresholds, any CBC inclusion lists, and limits
        // on number of CBCs to take from a DGE.
        List<DgeDescr> dgeDescrs = loadDataSetsPass1(dataSets);

        // Apply MIN_CELLS threshold to genes that have been seen.
        final GeneFiltererSorter geneFiltererSorter = new GeneFiltererSorter(MIN_CELLS);

        // update the gene-CBC bit matrix for each DGE to reflect the genes that have been filtered,
        // then remove any CBCs which have lost all their non-zero genes.
        filterGenesAndEmptyCells(dgeDescrs, geneFiltererSorter.filteredGeneIndices);

        final long numNonZeroElements = countNonZeroElements(dgeDescrs);
        writeDiscardedCellsFile(dgeDescrs);
        writeDgeHeader(dataSets);


        final List<String> cellBarcodes = new ArrayList<>();
        for (final DgeDescr dgeDescr : dgeDescrs) {
            cellBarcodes.addAll(dgeDescr.cellBarcodes);
        }
        final int[] cellSizes = new int[cellBarcodes.size()];

        // Write the integer DGE, if requested.  Even if integer DGE is not requested, it is still necessary
        // to read the input DGEs in order to determine the final cell sizes after thresholding, in order
        // to write the scaled DGE.
        final MatrixMarketWriter rawDgeWriter;
        if (RAW_DGE_OUTPUT_FILE != null) {
            rawDgeWriter = new MatrixMarketWriter(RAW_DGE_OUTPUT_FILE, MatrixMarketConstants.ElementType.integer,
                    geneFiltererSorter.getSortedGeneNames().size(), cellBarcodes.size(),
                    numNonZeroElements, geneFiltererSorter.getSortedGeneNames(), cellBarcodes,
                    MatrixMarketConstants.GENES, MatrixMarketConstants.CELL_BARCODES);
        } else {
            rawDgeWriter = null;
        }

        int cellIndexOffset = 0;
        int numFilteredElements = 0;
        for (final DgeDescr dgeDescr : dgeDescrs) {
            LOG.info("Counting cell sizes and writing raw DGE " + dgeDescr.file.getAbsolutePath());
            final RawLoadedDge dge = new RawLoadedDge(dgeDescr.file, geneEnumerator);
            dgeDescr.cellBarcodeIndexMap = new CellBarcodeIndexMap(dge.rawCellBarcode, dgeDescr.cellBarcodes, dgeDescr.name);
            for (final SparseDge.Triplet triplet: dge.rawTriplets) {
                final int geneIndex = geneFiltererSorter.getOutputGeneIndex(triplet.geneIndex);
                if (geneIndex < 0) {
                    // Gene was filtered by a threshold.
                    ++numFilteredElements;
                    continue;
                }
                if (dgeDescr.cellBarcodeIndexMap.oldToNewIndex[triplet.cellIndex] != -1) {
                    final int cellIndex = cellIndexOffset + dgeDescr.cellBarcodeIndexMap.oldToNewIndex[triplet.cellIndex];
                    cellSizes[cellIndex] += triplet.value;
                    if (rawDgeWriter != null) {
                        rawDgeWriter.writeTriplet(geneIndex, cellIndex, triplet.value);
                    }
                }
            }
            cellIndexOffset += dgeDescr.cellBarcodes.size();
        }
        try {
            if (rawDgeWriter != null) {
                rawDgeWriter.close();
            }
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }

        writeCellSizesFile(cellBarcodes, cellSizes);

        // Read the inputs again, and write scaled DGE, using cellSizes computed above.
        if (SCALED_DGE_OUTPUT_FILE != null) {
            final MatrixMarketWriter scaledDgeWriter =
                    new MatrixMarketWriter(SCALED_DGE_OUTPUT_FILE, MatrixMarketConstants.ElementType.real,
                            geneFiltererSorter.getSortedGeneNames().size(), cellBarcodes.size(),
                            numNonZeroElements, geneFiltererSorter.getSortedGeneNames(), cellBarcodes,
                            MatrixMarketConstants.GENES, MatrixMarketConstants.CELL_BARCODES);

            cellIndexOffset = 0;
            for (final DgeDescr dgeDescr : dgeDescrs) {
                LOG.info("Writing scaled DGE " + dgeDescr.file.getAbsolutePath());
                final RawLoadedDge dge = new RawLoadedDge(dgeDescr.file, geneEnumerator);
                for (final SparseDge.Triplet triplet: dge.rawTriplets) {
                    final int geneIndex = geneFiltererSorter.getOutputGeneIndex(triplet.geneIndex);
                    if (geneIndex < 0) {
                        // Gene was filtered by a threshold.
                        continue;
                    }
                    if (dgeDescr.cellBarcodeIndexMap.oldToNewIndex[triplet.cellIndex] != -1) {
                        final int cellIndex = cellIndexOffset + dgeDescr.cellBarcodeIndexMap.oldToNewIndex[triplet.cellIndex];
                        scaledDgeWriter.writeTriplet(geneIndex, cellIndex, triplet.value / (double) cellSizes[cellIndex]);
                    }
                }
                cellIndexOffset += dgeDescr.cellBarcodes.size();
            }
            try {
                scaledDgeWriter.close();
            } catch (IOException e) {
                throw new RuntimeIOException(e);
            }
        }

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
                        if (!fields[0].isEmpty())
                            ret.add(fields[0]);
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
        if (this.CELL_SIZE_OUTPUT_FILE!=null) IOUtil.assertFileIsWritable(this.CELL_SIZE_OUTPUT_FILE);
        if (this.RAW_DGE_OUTPUT_FILE!=null) IOUtil.assertFileIsWritable(this.RAW_DGE_OUTPUT_FILE);
        if (this.SCALED_DGE_OUTPUT_FILE!=null) IOUtil.assertFileIsWritable(this.SCALED_DGE_OUTPUT_FILE);
        if (this.DGE_HEADER_OUTPUT_FILE!=null) IOUtil.assertFileIsWritable(this.DGE_HEADER_OUTPUT_FILE);
        if (this.DISCARDED_CELLS_FILE!=null) IOUtil.assertFileIsWritable(this.DISCARDED_CELLS_FILE);

        final ArrayList<String> list = new ArrayList<>(1);

        if (RAW_DGE_OUTPUT_FILE == null && SCALED_DGE_OUTPUT_FILE == null) {
            list.add("At least one of RAW_DGE_OUTPUT_FILE and SCALED_DGE_OUTPUT_FILE should be set");
        }
        if (MIN_CELLS < 1) {
            list.add("MIN_CELLS must be >= 1.");
        }
        if (MIN_GENES < 1) {
            list.add("MIN_GENES must be >= 1.");
        }
        if (MIN_TRANSCRIPTS < 1) {
            list.add("MIN_TRANSCRIPTS must be >= 1.");
        }
        return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);
    }

    private class CellsPerGeneCounter {
        // gene index is from GeneEnumerator
        private  int[] cellsPerGene = null;
        void addCellCounts(final int[] countsToAdd) {
            if (cellsPerGene == null) {
                cellsPerGene = Arrays.copyOf(countsToAdd, countsToAdd.length);
            } else {
                if (cellsPerGene.length < countsToAdd.length) {
                    cellsPerGene = Arrays.copyOf(cellsPerGene, countsToAdd.length);
                }
                for (int i = 0; i < countsToAdd.length; ++i) {
                    cellsPerGene[i] += countsToAdd[i];
                }
            }
        }

        void addCellCounts(final SparseDge dge) {
            final int[] cellsPerGene = new int[geneEnumerator.getGenes().size()];
            for (final SparseDge.Triplet triplet : dge.getTriplets()) {
                ++cellsPerGene[triplet.geneIndex];
            }
            addCellCounts(cellsPerGene);
        }

    }

    private class GeneFiltererSorter {
        private final int[] geneIdMapping;
        private final List<String> sortedGeneNames;
        private final int[] filteredGeneIndices;

        public GeneFiltererSorter(final int minCellsPerGene) {
            final int[] cellsPerGene = cellsPerGeneCounter.cellsPerGene;
            geneIdMapping = new int[geneEnumerator.getGenes().size()];
            final Map<String, Integer> geneMap = new TreeMap<>();
            int numFilteredGenes = 0;
            for (int i = 0; i < geneIdMapping.length; ++i) {
                if (cellsPerGene[i] >= minCellsPerGene)
                    geneMap.put(geneEnumerator.getGeneName(i), i);
                else {
                    geneIdMapping[i] = -1;
                    ++numFilteredGenes;
                }
            }
            int numOutputGenes = 0;
            for (final Map.Entry<String, Integer> entry : geneMap.entrySet())
                geneIdMapping[entry.getValue()] = numOutputGenes++;
            sortedGeneNames = Collections.unmodifiableList(new ArrayList<>(geneMap.keySet()));
            filteredGeneIndices = new int[numFilteredGenes];
            int j = 0;
            for (int i = 0; i < geneIdMapping.length; ++i) {
                if (geneIdMapping[i] == -1) {
                    filteredGeneIndices[j++] = i;
                }
            }
        }

        public int getOutputGeneIndex(final int originalGeneIndex) {
            return geneIdMapping[originalGeneIndex];
        }

        public List<String> getSortedGeneNames() {
            return sortedGeneNames;
        }
    }

    private static class DgeBitSet {
        // One element for each cell barcode in a DGE
        // Ifthe ith bit is set, the ith gene (from GeneEnumerator) has non-zero expression
        final BitSet[] genesForCellBarcode;

        public DgeBitSet(final int numCellBarcodes, final int numGenes) {
            genesForCellBarcode = new BitSet[numCellBarcodes];
            for (int i = 0; i < numCellBarcodes; ++i) {
                genesForCellBarcode[i] = new BitSet(numGenes);
            }
        }
    }

    private void writeDgeHeader(final List<Map<String, Object>> dataSets) {
        if (DGE_HEADER_OUTPUT_FILE != null) {
            LOG.info("Writing " + DGE_HEADER_OUTPUT_FILE.getAbsolutePath());
            final List<File> inputDges = new ArrayList<>(dataSets.size());
            final List<String> prefixes = new ArrayList<>(dataSets.size());
            for (final Map<String, Object> dataSet : dataSets) {
                inputDges.add(new File((String)dataSet.get(YamlKeys.DatasetsKeys.PATH_KEY)));
                prefixes.add((String)getValueOrDefault(dataSet, YamlKeys.DatasetsKeys.NAME_KEY, ""));
            }
            new DgeHeaderCodec().encode(DGE_HEADER_OUTPUT_FILE, DgeHeaderMerger.mergeDgeHeaders(inputDges, prefixes, HEADER_STRINGENCY));
        }
    }

    private void writeCellSizesFile(final List<String> cellBarcodes, final int[] cellSizes) {
        if (cellBarcodes.size() != cellSizes.length) {
            throw new IllegalArgumentException("unpossible!");
        }
        if (CELL_SIZE_OUTPUT_FILE != null) {
            LOG.info("Writing " + CELL_SIZE_OUTPUT_FILE.getAbsolutePath());
            final CellSizeWriter cellSizeWriter = new CellSizeWriter(CELL_SIZE_OUTPUT_FILE);
            for (int i = 0; i < cellBarcodes.size(); ++i) {
                cellSizeWriter.writeSize(cellBarcodes.get(i), cellSizes[i]);
            }
            cellSizeWriter.close();
        }
    }

    private void writeDiscardedCellsFile(final List<DgeDescr> dgeDescrs) {
        if (DISCARDED_CELLS_FILE != null) {
            try {
                LOG.info("Writing " + DISCARDED_CELLS_FILE.getAbsolutePath());
                final BufferedWriter writer = IOUtil.openFileForBufferedWriting(DISCARDED_CELLS_FILE);
                for (final DgeDescr dgeDescr : dgeDescrs) {
                    for (final String cell_barcode: dgeDescr.discardedCells) {
                        writer.write(cell_barcode);
                        writer.newLine();
                    }
                }
                writer.close();
            } catch (IOException e) {
                throw new RuntimeIOException("Exception writing " + DISCARDED_CELLS_FILE.getAbsolutePath(), e);
            }
        }
    }

    private List<DgeDescr> loadDataSetsPass1(final List<Map<String, Object>> datasets) {
        if (datasets.size() == 0)
            throw new RuntimeException("List of datasets is empty in " + YAML.getAbsolutePath());

        // There can be the same prefix used multiple times with several DGE files.
        // In this case, we need to check that after filtering and prefixing any 2 DGEs associated with the same prefix contain non-overlapping prefixed cell barcode sets
        // We compare the prefixed cell barcode sets, since there can be a situation where one or more datasets for the same prefix do not have cell_barcode_path specified
        final ArrayList<DgeDescr> datasetDescrs = new ArrayList<>(datasets.size());

        for (final Map<String, Object> datasetMap : datasets) {
            final File dgePath = new File((String)getRequiredValue(datasetMap, YamlKeys.DatasetsKeys.PATH_KEY));
            final String cellBarcodePath = (String) getValueOrDefault(datasetMap, YamlKeys.DatasetsKeys.CELL_BARCODE_PATH_KEY, null);
            String prefix = (String)getValueOrDefault(datasetMap, YamlKeys.DatasetsKeys.NAME_KEY, "");
            final SparseDge dge = new SparseDge(dgePath, geneEnumerator);
            LOG.info(String.format("Loaded %d cells from %s", dge.getNumCells(), dgePath.getAbsolutePath()));

            // If specified, this file contains actual un-prefixed cell barcodes for the cells to be retained from this DGE file
            if (cellBarcodePath != null) {
                final File cellBarcodeFile = new File(cellBarcodePath);
                Set<String> unprefixedSelectedCells = loadSelectedCellsLists(Collections.singletonList(cellBarcodeFile));
                dge.retainOnlyTheseCells(unprefixedSelectedCells);
                LOG.info(String.format("After applying %s, %d cells remain.", cellBarcodeFile.getAbsolutePath(), dge.getNumCells()));
            }

            if (!prefix.isEmpty())
                dge.prefixCellBarcodes(prefix + "_");
            Integer cell_count = (Integer)getValueOrDefault(datasetMap, YamlKeys.DatasetsKeys.CELL_COUNT_KEY, 0);
            if (cell_count > 0 && cell_count < dge.getNumCells())
                dge.discardSmallestCells(cell_count);

            if (selectedCells != null) {
                int numCellsBefore = dge.getNumCells();
                dge.retainOnlyTheseCells(selectedCells);
                LOG.info(String.format("%d cells removed after applying global CELL_BC_FILES.", numCellsBefore - dge.getNumCells()));
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

            cellsPerGeneCounter.addCellCounts(dge);
            final List<String> cellBarcodes = dge.getCellBarcodes();
            // Check that this cell barcode set does not overlap any other one associated with the same prefix
            datasetDescrs.stream().filter(descr -> prefix.equals(descr.name)).collect(Collectors.toList()).forEach(descr -> {
                List<String> commonCellBarcodes = descr.cellBarcodes.stream()
                    .filter(cellBarcodes::contains)
                    .collect(Collectors.toList());
                if (!commonCellBarcodes.isEmpty()) {
                    throw new RuntimeException(String.format("DGE files %s and %s both have the prefix %s, and yet they contain %d cell barcodes in common", descr.file.getAbsolutePath(), dge.getFile().getAbsolutePath(), prefix, commonCellBarcodes.size()));
                }
            });
            final DgeBitSet dgeBitSet = new DgeBitSet(cellBarcodes.size(), geneEnumerator.getNumGenes());
            for (final SparseDge.Triplet triplet : dge.getTriplets()) {
                if (triplet.value == 0) {
                    throw new IllegalStateException("Unpossible!");
                }
                dgeBitSet.genesForCellBarcode[triplet.cellIndex].set(triplet.geneIndex);
            }
            datasetDescrs.add(new DgeDescr(dgePath, prefix, cellBarcodes, dge.getDiscardedCells(), dgeBitSet));
        }

        return datasetDescrs;
    }

    private void filterGenesAndEmptyCells(final List<DgeDescr> dgeDescrs, final int[] filteredGeneIndices) {
        for (final DgeDescr dgeDescr : dgeDescrs){
            // Clear all the filtered genes from the bit set
            final DgeBitSet dgeBitSet = dgeDescr.dgeBitSet;
            for (final BitSet bitSet : dgeBitSet.genesForCellBarcode) {
                for (final int geneIndex : filteredGeneIndices) {
                    if (geneIndex < bitSet.length()) {
                        bitSet.clear(geneIndex);
                    }
                }
            }
            final ArrayList<String> cellBarcodes = dgeDescr.cellBarcodes;
            if (dgeBitSet.genesForCellBarcode.length != dgeDescr.cellBarcodes.size()) {
                throw new IllegalStateException("Unpossible");
            }
            // Find any cell barcodes that now have zero expression
            // Start from the end of the array so that elements can be removed without messing up indices
            // of previous cell barcodes.
            for (int i = cellBarcodes.size() - 1; i >= 0; --i) {
                if (dgeBitSet.genesForCellBarcode[i].isEmpty()) {
                    dgeDescr.discardedCells.add(cellBarcodes.get(i));
                    cellBarcodes.remove(i);
                }
            }
        }
    }

    private long countNonZeroElements(final List<DgeDescr> dgeDescrs) {
        long numNonZeroElements = 0L;
        for (final DgeDescr dgeDescr: dgeDescrs) {
            numNonZeroElements += Arrays.stream(dgeDescr.dgeBitSet.genesForCellBarcode).mapToLong(BitSet::cardinality).sum();
        }
        return numNonZeroElements;
    }


    private Object getRequiredValue(final Map<String, Object> map, final String key) {
        final Object ret = map.get(key);
        if (ret == null)
            throw new RuntimeException(YAML.getAbsolutePath() + " does not contain key " + key);
        return ret;
    }

    private Object getValueOrDefault(final Map<String, Object> map, final String key, final Object defaultValue) {
        final Object ret = map.get(key);
        if (ret == null)
            return defaultValue;
        else
            return ret;
    }

    private static final class CellBarcodeIndexMap {
        final int[] oldToNewIndex;

        public CellBarcodeIndexMap(final String[] rawCellBarcodes, final ArrayList<String> finalCellBarcodes, final String name) {
            oldToNewIndex = new int[rawCellBarcodes.length];
            final String prefix = (name == null || name.isEmpty()? "" : name + "_");
            for (int i = 0; i < rawCellBarcodes.length; ++i) {
                oldToNewIndex[i] = finalCellBarcodes.indexOf(prefix + rawCellBarcodes[i]);
            }
        }
    }
}
