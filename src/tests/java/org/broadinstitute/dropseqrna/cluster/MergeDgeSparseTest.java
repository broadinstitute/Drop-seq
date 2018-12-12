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
import htsjdk.samtools.util.TestUtil;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeaderMerger;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.tools.DGEMatrix;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import org.yaml.snakeyaml.Yaml;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

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
        final File tempDir = TestUtil.getTempDirectory("MergeDgeSparseTest.", ".tmp");
        tempDir.deleteOnExit();

        final Yaml yamlConverter = new Yaml();
        final Map yamlMap = (Map)yamlConverter.load(IOUtil.openFileForReading(YAML));
        final List datasets = (List)yamlMap.get(MergeDgeSparse.YamlKeys.DATASETS_KEY);

        if (cell_count != null) {
            for (int i = 0; i < cell_count.length; ++i) {
                ((Map)datasets.get(i)).put(MergeDgeSparse.YamlKeys.DatasetsKeys.CELL_COUNT_KEY, cell_count[i]);
            }
        }
        if (prefix != null) {
            for (int i = 0; i < prefix.length; ++i) {
                ((Map)datasets.get(i)).put(MergeDgeSparse.YamlKeys.DatasetsKeys.NAME_KEY, prefix[i]);
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
    }

    @DataProvider(name="testBasicDataProvider")
    private Object[][] testBasicDataProvider() throws IOException {
        final List<File> selectedCellsFiles = new ArrayList<File>();
        for (final Path path: Files.newDirectoryStream(TEST_DATA_DIR.toPath(), "selected_cells.*.txt*")) {
            selectedCellsFiles.add(path.toFile());
        }
        return new Object[][]{
                {"basic", 109, 297, null, new String[]{"A", "B", "C"}, null, null, 0, null, DgeHeaderMerger.Stringency.STRICT, null},
                {"cellCount", 109, 248, new int[]{50, 0, 200}, new String[]{"A", "B", "C"}, null, null, 0, null, DgeHeaderMerger.Stringency.STRICT, null},
                {"noPrefix", 109, 297, null, null, null, null, 0, null, DgeHeaderMerger.Stringency.STRICT, null},
                {"filterGenes", 108, 297, null, new String[]{"A", "B", "C"}, new String[]{"^mt-"}, null, 0, null, DgeHeaderMerger.Stringency.STRICT, null},
                {"minCells", 74, 297, null, new String[]{"A", "B", "C"}, null, 2, 0, null, DgeHeaderMerger.Stringency.STRICT, null},
                {"minGenes", 109, 171, null, new String[]{"A", "B", "C"}, null, null, 10, null, DgeHeaderMerger.Stringency.STRICT, null},
                {"minTranscripts", 109, 236, null, new String[]{"A", "B", "C"}, null, null, 0, 10, DgeHeaderMerger.Stringency.STRICT, null},
                {"selectedCells", 109, 9, null, new String[]{"A", "B", "C"}, null, null, 0, null, DgeHeaderMerger.Stringency.STRICT, selectedCellsFiles},
        };
    }
}
