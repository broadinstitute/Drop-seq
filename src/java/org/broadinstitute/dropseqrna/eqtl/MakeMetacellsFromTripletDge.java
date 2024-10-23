/*
 * MIT License
 *
 * Copyright 2024 Broad Institute
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
package org.broadinstitute.dropseqrna.eqtl;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.matrixmarket.MatrixMarketConstants;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.TabbedInputParser;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.BufferedReader;
import java.io.File;
import java.io.PrintStream;
import java.util.*;

@CommandLineProgramProperties(
        summary = """
        Create a metacell file from a triplet DGE (matrix.mtz.gz, features.tsv.gz, barcodes.tsv.gz) and
        a tabular mapping file that details which cell barcodes belong to which metacell.
        """,
        oneLineSummary = "Create a metacell file from a triplet DGE and a tabular mapping file.",
        programGroup = DropSeq.class
)
public class MakeMetacellsFromTripletDge
extends CommandLineProgram {

    @Argument(doc = "The triplet DGE matrix file (typically matrix.mtx.gz).", shortName = "M")
    public File MATRIX;

    @Argument(doc = "The triplet DGE features file.  If not specified, a file named features.tsv.gz in the same " +
            "directory as MATRIX will be used", shortName = "F", optional = true)
    public File FEATURES;

    @Argument(doc = "The triplet DGE barcodes file.  If not specified, a file named barcodes.tsv.gz in the same " +
            "directory as MATRIX will be used", shortName = "B", optional = true)
    public File BARCODES;

    @Argument(doc = "The tabular mapping file with header that details which cell barcodes belong to which metacell.")
    public File MAPPING;

    @Argument(doc = "The output metacell file.  If ends with .gz, output will bw compressed.",
            shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME)
    public File OUTPUT;

    @Argument(doc = "The columns in the MAPPING file that should be used to group the cell barcodes into metacells.",
            minElements = 1, shortName = "G")
    public List<String> GROUP_COLUMNS;

    @Argument(doc="The name of the column in the mapping file that contains the cell barcode.")
    public String CELL_BARCODE_COLUMN = "CELL_BARCODE";
    @Argument(doc="The name of the column in the mapping file that contains the prefix to be prepended to " +
            "the CBC to match the values in BARCODES file.  Set to null to disable prefixing.", optional = true)
    public String PREFIX_COLUMN = "PREFIX";

    @Argument(doc="This string is used to join the prefix to the CBC in the MAPPING file")
    public String PREFIX_SEPARATOR = "_";

    @Argument(doc="This string is used to join the values in the group columns to create the metacell name.")
    public String GROUP_SEPARATOR = ":";

    @Argument( doc = "Value to use for missing group values.  Default: skip rows in MAPPING file that have an empty " +
            "value for one or more of the GROUP_COLUMNS.",optional = true)
    public String MISSING_GROUP_VALUE;

    public static final String DEFAULT_BARCODES_FILENAME = "barcodes.tsv.gz";
    public static final String DEFAULT_FEATURES_FILENAME = "features.tsv.gz";

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(MATRIX);
        IOUtil.assertFileIsReadable(MAPPING);
        IOUtil.assertFileIsReadable(BARCODES);
        IOUtil.assertFileIsReadable(FEATURES);
        IOUtil.assertFileIsWritable(OUTPUT);
        final List<String> barcodes = readList(BARCODES);
        final List<String> features = readList(FEATURES);
        final ParsedMapping parsedMapping = createBarcodeToMetacellMap(barcodes);

        // Open the matrix market file, do some validation and advance to the triplets
        final BufferedReader reader = IOUtil.openFileForBufferedReading(MATRIX);
        String line = readLine(reader, false);
        if (!MatrixMarketConstants.MM_HEADER_INT.equals(line)) {
            throw new IllegalArgumentException(String.format("Matrix file %s does not start with expected header line: %s\nactual: %s",
                    MATRIX, MatrixMarketConstants.MM_HEADER_INT, line));
        }
        // Skip the header comments
        while ((line = readLine(reader, false)).startsWith(MatrixMarketConstants.MM_COMMENT_LINE_START)) {
            // Skip
        }
        final long[] triplet = new long[3];
        // Read the dimensions and validate
        parseTriplet(line, triplet);
        if (triplet[0] != features.size() || triplet[1] != barcodes.size()) {
            throw new IllegalArgumentException(String.format("Matrix file %s has dimensions %d x %d, expected %d x %d based on %s and %s",
                    MATRIX, triplet[0], triplet[1], features.size(), barcodes.size(), FEATURES, BARCODES));
        }
        // Accumulate the counts for each metacell a gene at a time
        final long[] expressionForGene = new long[parsedMapping.metacellNames.length];
        // read the first triplet and prepare to read the rest
        readTriplet(reader, triplet, false);
        long geneIndex = triplet[0];
        int metacellIndex = parsedMapping.cellBarcodeMetacellIndex[(int)triplet[1]];
        if (metacellIndex != -1) {
            expressionForGene[metacellIndex] += triplet[2];
        }

        PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(OUTPUT));
        // write header
        out.println("#"+getCommandLine());
        out.println(CreateMetaCells.GENE_HEADER + "\t" + StringUtil.join("\t", parsedMapping.metacellNames));
        // Read all the triplets for the same gene, and tally the counts for each metacell.
        while (readTriplet(reader, triplet, true)) {
            if (triplet[0] != geneIndex) {
                if (triplet[0] < geneIndex) {
                    throw new IllegalArgumentException(
                            String.format("Gene indices not in order i %s: %d < %d", MATRIX, triplet[0], geneIndex));
                }
                // Write out the tallied counts for the gene
                writeGene(out, features.get((int)geneIndex-1), expressionForGene);
                // Reset the tallies for the next gene
                Arrays.fill(expressionForGene, 0);
                geneIndex = triplet[0];
            }
            metacellIndex = parsedMapping.cellBarcodeMetacellIndex[(int)triplet[1]];
            if (metacellIndex != -1) {
                expressionForGene[metacellIndex] += triplet[2];
            }        }
        // Write out the tallied counts for the last gene
        writeGene(out, features.get((int)geneIndex-1), expressionForGene);
        out.close();
        CloserUtil.close(reader);

        return 0;
    }

    @Override
    protected String[] customCommandLineValidation() {
        final ArrayList<String> list = new ArrayList<>(1);
        if (GROUP_COLUMNS.isEmpty()) {
            list.add("GROUP_COLUMNS must contain at least one column name.");
        }
        if (BARCODES == null) {
            BARCODES = new File(MATRIX.getParent(), DEFAULT_BARCODES_FILENAME);
        }
        if (FEATURES == null) {
            FEATURES = new File(MATRIX.getParent(), DEFAULT_FEATURES_FILENAME);
        }
        return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);

    }

    // For both barcodes and features, we only care about first column, so the fact there are multiple
    // columns in features file is not a problem.
    private List<String> readList(final File inputFile) {
        final TabbedInputParser parser = new TabbedInputParser(false, inputFile);
        final List<String> ret = new ArrayList<>();
        for (final String[] row: parser) {
            ret.add(row[0]);
        }
        CloserUtil.close(parser);
        return ret;
    }

    private record ParsedMapping(String[] metacellNames, int[] cellBarcodeMetacellIndex) {
    }

    private ParsedMapping createBarcodeToMetacellMap(final List<String> barcodes) {
        final int[] cellBarcodeMetacellIndices = new int[barcodes.size()+1]; // column numbers in MM file are 1-based
        Arrays.fill(cellBarcodeMetacellIndices, -1);
        // Use a LinkedHashMap to preserve the order of first appearance of a metacell in the mapping file.
        final Map<String, Integer> metacellToColumnMap = new LinkedHashMap<>();
        // Read the mapping file and with a TabbedTextFileWithHeaderParser, determine the CBC using CELL_BARCODE_COLUMN
        // and PREFIX column, and the metacell name using GROUP_COLUMNS.
        final TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(MAPPING);
        final Set<String> columnNames = parser.getColumnNames();
        if (!columnNames.contains(CELL_BARCODE_COLUMN)) {
            throw new IllegalArgumentException("Mapping file does not contain column " + CELL_BARCODE_COLUMN);
        }
        if (PREFIX_COLUMN != null && !columnNames.contains(PREFIX_COLUMN)) {
            throw new IllegalArgumentException("Mapping file does not contain column " + PREFIX_COLUMN);
        }
        for (final String groupColumn: GROUP_COLUMNS) {
            if (!columnNames.contains(groupColumn)) {
                throw new IllegalArgumentException("Mapping file does not contain column " + groupColumn);
            }
        }
        for (final TabbedTextFileWithHeaderParser.Row row: parser) {
            final String cellBarcode = getPrefixedBarcode(row);
            int barcodeIndex = barcodes.indexOf(cellBarcode);
            if (barcodeIndex == -1) {
                // Should this be an error?
                continue;
            }
            if (cellBarcodeMetacellIndices[barcodeIndex+1] != -1) {
                throw new IllegalArgumentException("Duplicate cell barcode in mapping file: " + cellBarcode);
            }
            List<String> group_values = GROUP_COLUMNS.stream().map(row::getField).toList();
            if (MISSING_GROUP_VALUE == null) {
                if (group_values.contains(null)) {
                    continue; // skip this row
                }
            } else {
                group_values = group_values.stream().map(v -> v == null ? MISSING_GROUP_VALUE : v).toList();
            }
            final String metacell = StringUtil.join(GROUP_SEPARATOR, group_values);
            final Integer metacellIndex = metacellToColumnMap.computeIfAbsent(metacell, k -> metacellToColumnMap.size());
            cellBarcodeMetacellIndices[barcodeIndex+1] = metacellIndex;
        }
        CloserUtil.close(parser);
        return new ParsedMapping(metacellToColumnMap.keySet().toArray(new String[0]), cellBarcodeMetacellIndices);
    }

    private String getPrefixedBarcode(final TabbedTextFileWithHeaderParser.Row row) {
        final String cellBarcode = row.getField(CELL_BARCODE_COLUMN);
        if (PREFIX_COLUMN == null) {
            return cellBarcode;
        }
        return row.getField(PREFIX_COLUMN) + PREFIX_SEPARATOR + cellBarcode;
    }

    private String readLine(final BufferedReader reader, final boolean eofOK) {
        try {
            final String line = reader.readLine();
            if (line == null && !eofOK) {
                throw new IllegalArgumentException("Unexpected EOF in matrix file " + MATRIX);
            }
            return line;
        } catch (Exception e) {
            throw new RuntimeIOException("Error reading matrix file " + MATRIX, e);
        }
    }

    private boolean readTriplet(final BufferedReader reader, final long[] triplet, final boolean eofOK) {
        final String line = readLine(reader, eofOK);
        if (line == null) {
            return false;
        }
        parseTriplet(line, triplet);
        return true;
    }

    private void parseTriplet(final String line, final long[] triplet) {
        final String[] parts = StringUtils.split(line);
        if (parts.length != 3) {
            throw new IllegalArgumentException("Expected 3 parts in triplet line, got " + parts.length);
        }
        triplet[0] = Long.parseLong(parts[0]);
        triplet[1] = Long.parseLong(parts[1]);
        triplet[2] = Long.parseLong(parts[2]);
        if (triplet[0] < 1 || triplet[1] < 1 || triplet[2] < 0) {
            throw new IllegalArgumentException("Invalid triplet: " + line);
        }
    }

    private void writeGene(final PrintStream out, final String gene, final long[] counts) {
        if (Arrays.stream(counts).allMatch(c -> c == 0)) {
            // Skip genes with no expression for selected CBCs
            return;
        }
        out.print(gene);
        for (long count: counts) {
            out.print("\t");
            out.print(count);
        }
        out.println();
    }

    public static void main(final String[] args) {
        new MakeMetacellsFromTripletDge().instanceMainWithExit(args);
    }
}
