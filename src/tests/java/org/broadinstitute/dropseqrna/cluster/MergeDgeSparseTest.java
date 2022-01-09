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
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeaderMerger;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.tools.DGEMatrix;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.tools.MatrixTransformFactory;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.yaml.snakeyaml.Yaml;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

public class MergeDgeSparseTest {
    private static final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/cluster");
    private static final File YAML = new File(TEST_DATA_DIR, "test.yaml");

    @SuppressWarnings("unchecked")
    @Test(dataProvider = "testBasicDataProvider")
    public void testBasic(final String testName,
                          final int expectedNumGenes,
                          final int expectedNumCells,
                          final int[] cell_count,
                          final String[] prefix,
                          final String[] filteredGeneRe,
                          final Integer min_cells,
                          final Integer min_genes,
                          final Integer min_transcripts,
                          final DgeHeaderMerger.Stringency headerStringency,
                          final List<File> selectedCellsFiles) throws IOException {
        final File tempDir = Files.createTempDirectory("MergeDgeSparseTest.").toFile();
        tempDir.deleteOnExit();

        final Yaml yamlConverter = new Yaml();
        final Map<String, Object> yamlMap = (Map<String, Object>)yamlConverter.load(IOUtil.openFileForReading(YAML));
        final List<Map<String, Object>> datasets = (List<Map<String, Object>>)yamlMap.get(MergeDgeSparse.YamlKeys.DATASETS_KEY);

        if (cell_count != null) {
            for (int i = 0; i < cell_count.length; ++i) {
                (datasets.get(i)).put(MergeDgeSparse.YamlKeys.DatasetsKeys.CELL_COUNT_KEY, cell_count[i]);
            }
        }
        if (prefix != null) {
            for (int i = 0; i < prefix.length; ++i) {
                (datasets.get(i)).put(MergeDgeSparse.YamlKeys.DatasetsKeys.NAME_KEY, prefix[i]);
            }
        }
        final File outputYaml = new File(tempDir, "test.yaml");
        outputYaml.deleteOnExit();
        final BufferedWriter outputYamlWriter = IOUtil.openFileForBufferedWriting(outputYaml);
        yamlConverter.dump(yamlMap, outputYamlWriter);
        outputYamlWriter.close();

        final File cellSizeOutputFile = new File(tempDir, "test.cell_size.txt");    cellSizeOutputFile.deleteOnExit();
        final File dgeHeaderOutputFile = new File(tempDir, "test.dge_header.txt");  dgeHeaderOutputFile.deleteOnExit();
        final File rawDgeOutputFile = new File(tempDir, "test.raw.dge.txt");        rawDgeOutputFile.deleteOnExit();
        final File scaledDgeOutputFile = new File(tempDir, "test.scaled.dge.txt");  scaledDgeOutputFile.deleteOnExit();
        final File discardedCellsOutputFile = new File(tempDir, "test.discarded_cells.txt"); discardedCellsOutputFile.deleteOnExit();

        final MergeDgeSparse merger = new MergeDgeSparse();
        merger.YAML = outputYaml;
        merger.CELL_SIZE_OUTPUT_FILE = cellSizeOutputFile;
        if (prefix != null) {
            // If not setting prefixes, can't create merged header
            merger.DGE_HEADER_OUTPUT_FILE = dgeHeaderOutputFile;
        }
        merger.RAW_DGE_OUTPUT_FILE = rawDgeOutputFile;
        merger.SCALED_DGE_OUTPUT_FILE = scaledDgeOutputFile;
        merger.HEADER_STRINGENCY = headerStringency;

        if (min_cells != null) {
            merger.MIN_CELLS = min_cells;
        }
        if (min_genes != null) {
            merger.MIN_GENES = min_genes;
        }
        if (min_transcripts != null) {
            merger.MIN_TRANSCRIPTS = min_transcripts;
        }
        if (filteredGeneRe != null) {
            merger.FILTERED_GENE_RE = Arrays.asList(filteredGeneRe);
        } else {
            merger.FILTERED_GENE_RE = Collections.emptyList();
        }
        merger.CELL_BC_FILE = selectedCellsFiles;
        merger.DISCARDED_CELLS_FILE = discardedCellsOutputFile;

        Assert.assertEquals(merger.doWork(), 0);

        final DGEMatrix rawMatrix = DGEMatrix.parseFile(rawDgeOutputFile);
        final DGEMatrix scaledMatrix = DGEMatrix.parseFile(scaledDgeOutputFile);
        Assert.assertEquals(rawMatrix.getGenes().size(), scaledMatrix.getGenes().size(), "nrow(raw matrix) != nrow(scaled matrix)");
        Assert.assertEquals(rawMatrix.getCellBarcodes().size(), scaledMatrix.getCellBarcodes().size(), "ncol(raw matrix) != ncol(scaled matrix)");
        Assert.assertEquals(rawMatrix.getGenes().size(), expectedNumGenes);
        Assert.assertEquals(rawMatrix.getCellBarcodes().size(), expectedNumCells);
        rawMatrix.applyTransform(MatrixTransformFactory.normalizeColumns());
        final double[][] rawExpressionMatrix = rawMatrix.getExpressionMatrix();
        final double[][] scaledExpressionMatrix = scaledMatrix.getExpressionMatrix();
        for (int i = 0; i < rawExpressionMatrix.length; ++i) {
            for (int j = 0; j < rawExpressionMatrix[i].length; ++j) {
                Assert.assertEquals(scaledExpressionMatrix[i][j], rawExpressionMatrix[i][j], 0.000001);
            }
        }
    }

    @DataProvider(name="testBasicDataProvider")
    private Object[][] testBasicDataProvider() throws IOException {
        final List<File> selectedCellsFiles = new ArrayList<>();
        for (final Path path: Files.newDirectoryStream(TEST_DATA_DIR.toPath(), "selected_cells.*.txt*")) {
            selectedCellsFiles.add(path.toFile());
        }
        return new Object[][]{
                {"basic", 85, 297, null, new String[]{"A", "B", "C"}, null, null, 1, null, DgeHeaderMerger.Stringency.STRICT, null},
                {"cellCount", 82, 248, new int[]{50, 0, 200}, new String[]{"A", "B", "C"}, null, null, 1, null, DgeHeaderMerger.Stringency.STRICT, null},
                {"noPrefix", 85, 297, null, null, null, null, 1, null, DgeHeaderMerger.Stringency.STRICT, null},
                {"filterGenes", 85, 297, null, new String[]{"A", "B", "C"}, new String[]{"^mt-"}, null, 1, null, DgeHeaderMerger.Stringency.STRICT, null},
                {"minCells", 74, 297, null, new String[]{"A", "B", "C"}, null, 2, 1, null, DgeHeaderMerger.Stringency.STRICT, null},
                {"minGenes", 84, 171, null, new String[]{"A", "B", "C"}, null, null, 10, null, DgeHeaderMerger.Stringency.STRICT, null},
                {"minTranscripts", 84, 236, null, new String[]{"A", "B", "C"}, null, null, 1, 10, DgeHeaderMerger.Stringency.STRICT, null},
                {"selectedCells", 47, 9, null, new String[]{"A", "B", "C"}, null, null, 1, null, DgeHeaderMerger.Stringency.STRICT, selectedCellsFiles},
        };
    }

    /**
     * When a gene is filtered because of MIN_CELLS, it can cause a cell barcode's expression to go to 0.  If that
     * happens, that cell barcode should not appear in the output.
     * The cells are sorted in descending size order, so the test parameters cause the cell to be discarded (cell 2)
     * to have the largest, smallest, or middle size.
     * @param c1Expression expression for cell 1 in genes 1 and 3
     * @param c2Expression expression for cell 2 in gene 2
     * @param c3Expression expression for cell 3 in genes 1 and 3
     */
    @Test(dataProvider = "testMinCellsCellDroppingProvider")
    public void testMinCellsCellDropping(final String testName, final int c1Expression, final int c2Expression, final int c3Expression)
            throws IOException {
        final String cellBarcodeToBeFiltered = "CCCC";
        final String prefix = "prefix";
        final String prefixedCellBarcodeToBeFiltered = prefix + "_" + cellBarcodeToBeFiltered;
        final File denseDgeFile = TestUtils.getTempReportFile("MergeDgeSparseTest.", ".digital_expression.txt");
        PrintStream denseDgePrintStream = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(denseDgeFile));
        denseDgePrintStream.println(StringUtil.join("\t", "GENE", "AAAA", cellBarcodeToBeFiltered, "GGGG"));
        final int numGenes = 3;
        for (int i = 0; i < numGenes; ++i) {
            final int thisGeneC1Expression;
            final int thisGeneC2Expression;
            final int thisGeneC3Expression;
            if (i == 1) {
                thisGeneC1Expression = 0;
                thisGeneC3Expression = 0;
                thisGeneC2Expression = c2Expression;
            } else {
                thisGeneC1Expression = c1Expression;
                thisGeneC3Expression = c3Expression;
                thisGeneC2Expression = 0;
            }
            denseDgePrintStream.println(StringUtil.join("\t", "GENE" + (i+1),
                    thisGeneC1Expression, thisGeneC2Expression, thisGeneC3Expression));
        }
        denseDgePrintStream.close();

        final Map<String, Object> yamlMap = new HashMap<>();
        final Map<String, Object> datasetMap = new HashMap<>();
        yamlMap.put(MergeDgeSparse.YamlKeys.DATASETS_KEY, Collections.singletonList(datasetMap));
        datasetMap.put(MergeDgeSparse.YamlKeys.DatasetsKeys.NAME_KEY, prefix);
        datasetMap.put(MergeDgeSparse.YamlKeys.DatasetsKeys.PATH_KEY, denseDgeFile.getPath());
        final File outputYaml = TestUtils.getTempReportFile("MergeDgeSparseTest.", ".yaml");
        final BufferedWriter outputYamlWriter = IOUtil.openFileForBufferedWriting(outputYaml);
        new Yaml().dump(yamlMap, outputYamlWriter);
        outputYamlWriter.close();

        final MergeDgeSparse merger = new MergeDgeSparse();
        merger.YAML = outputYaml;
        merger.RAW_DGE_OUTPUT_FILE = TestUtils.getTempReportFile("MergeDgeSparseTest.", ".mtx");
        merger.MIN_CELLS = 2;
        merger.MIN_GENES = 1;
        merger.MIN_TRANSCRIPTS = 1;
        merger.FILTERED_GENE_RE = Collections.emptyList();
        merger.DISCARDED_CELLS_FILE = TestUtils.getTempReportFile("MergeDgeSparseTest.", ".discarded_cells.txt");
        Assert.assertEquals(merger.doWork(), 0);

        final DGEMatrix rawMatrix = DGEMatrix.parseFile(merger.RAW_DGE_OUTPUT_FILE);
        Assert.assertEquals(rawMatrix.getCellBarcodes().size(), 2);
        Assert.assertEquals(rawMatrix.getGenes().size(), 2);
        Assert.assertFalse(rawMatrix.getCellBarcodes().contains(prefixedCellBarcodeToBeFiltered));
        Assert.assertFalse(rawMatrix.getGenes().contains("GENE2"));
        final List<String> discardedCells = StreamSupport.stream(
                IOUtil.readLines(merger.DISCARDED_CELLS_FILE).spliterator(), false).collect(Collectors.toList());
        Assert.assertEquals(discardedCells.size(), 1);
        Assert.assertEquals(discardedCells.get(0), prefixedCellBarcodeToBeFiltered);
    }
    @DataProvider(name="testMinCellsCellDroppingProvider")
    private Object[][] testMinCellsCellDroppingProvider() {
        return new Object[][]{
                {"largest", 10, 100, 9},
                {"smallest", 100, 10, 200},
                {"middle", 2, 10, 200},
        };
    }
}
