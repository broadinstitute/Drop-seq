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
package org.broadinstitute.dropseqrna.barnyard.digitalexpression;

import htsjdk.samtools.Defaults;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.matrixmarket.MatrixMarketConstants;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintWriter;
import org.yaml.snakeyaml.Yaml;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.TabbedInputParser;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.StreamSupport;

/**
 * Read a YAML manifest containing a list of DGEs and some configuration parameters.
 * Produce a Matrix Market sparse DGE, with {1-based row number, 1-based col number, value} tab-separated triplets
 * sorted by {row, column}.
 * Each DGE will be opened in order to determine the CBCs to retain, possibly filtering by a per-DGE
 * allow list.  This is enough information to write barcodes.tsv.gz file, and to determine the output
 * indices of the CBCs.
 * The DGEs will then be pushed into a PriorityQueue that will sort based on
 * {gene, ordinal of DGE in manifest} and write the triplets to a temporary, Snappy-compressed file.
 * The gene names will be learned as part of this process, and written to genes.tsv.gz.
 * The total number of non-zero matrix entries will be counted.
 * The final output will be written with the correct header and the count of rows, columns, and non-zero entries,
 * followed by the triplets read from the temporary file.
 *
 * Note that it is possible with CBC filtering to have a gene with 0 expression in the output matrix.
 */
@CommandLineProgramProperties(
        summary = """
                Read a YAML containing a list of DGEs and some configuration parameters,\s
                and produce a Matrix Market sparse DGE, with triplets sorted by {row, column}.\n
                Note that it is possible with CBC filtering to have a gene with 0 expression in the output matrix.
                """,
        oneLineSummary = "Merge DGEs into a sparse Matrix Market DGE",
        programGroup = DropSeq.class
)
public class MakeTripletDge
extends AbstractTripletDgeWriterClp {
    @Argument(shortName = "M", doc= """
            yaml input file containing list of DGEs to be merged.
            The file is expected to contain a 'datasets' list.  Each element of the list will contain:
            
            dge: the location of the tabular DGE, which must be lexically sorted by gene. (required)
            prefix: a prefix to prepend to each CBC.  (default: no prefix)
            barcode_list: File containing non-prefixed cell barcodes to include.  (default: include all CBCs in DGE)
            
            The DGEs will be merged in the order they appear in the list.
            CBCs (with optional prefix) must be unique.
            """)
    public File MANIFEST;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc="Output file in Matrix Market format.  Typically this is matrix.mtx.gz")
    public File OUTPUT;

    @Argument(doc="For a DGE with prefix specified in manifest, use this string to join the prefix to the CBC")
    public String PREFIX_SEPARATOR = "_";

    @Argument(doc="Dense matrix values that match this string will be considered empty.  Note that the DGEs are" +
            "not parsed as numbers, just treated as tab-separated strings.")
    public String ZERO = "0";

    @Argument(doc="Written to the header of the output file.")
    public MatrixMarketConstants.ElementType ELEMENT_TYPE = MatrixMarketConstants.ElementType.integer;

    public static void main(final String[] args) {
        new MakeTripletDge().instanceMainWithExit(args);
    }

    @Override
    protected int doWork() {
        final Map<String, Object> manifest = loadManifest(MANIFEST);
        final List<ParsedYamlDge> parsedDges = parseYamlDges(manifest, MANIFEST);
        final List<TabularDgeStream> dgeStreams = new ArrayList<>(parsedDges.size());
        final List<String> outputBarcodes = new ArrayList<>();
        for (int i = 0; i < parsedDges.size(); ++i) {
            final ParsedYamlDge parsedYamlDge = parsedDges.get(i);
            final TabularDgeStream tabularDgeStream = new TabularDgeStream(parsedYamlDge.dge, parsedYamlDge.prefix,
                    parsedYamlDge.barcodeList, parsedYamlDge.barcodeColumn, i, outputBarcodes.size());
            if (!tabularDgeStream.hasNext()) {
                throw new IllegalArgumentException("DGE " + parsedYamlDge.dge + " is empty");
            }
            dgeStreams.add(tabularDgeStream);
            outputBarcodes.addAll(tabularDgeStream.prefixedBarcodes);
        }
        final long numUniqueBarcodes = outputBarcodes.stream().distinct().count();
        if (numUniqueBarcodes != outputBarcodes.size()) {
            throw new IllegalArgumentException("Output CBCs are not unique");
        }
        if (outputBarcodes.isEmpty()) {
            throw new IllegalArgumentException("No CBCs to output");
        }
        dgeStreams.forEach(TabularDgeStream::clearBarcodes);
        final PriorityQueue<TabularDgeStream> queue = new PriorityQueue<>(dgeStreams.size(), new TabularDgeStreamComparator());
        queue.addAll(dgeStreams);
        final List<String> genes = new ArrayList<>();
        String currGene = "";  // Empty string sorts before everything and avoids null check
        long nonZeroCount = 0;
        final TempStreamFactory tempStreamFactory = new TempStreamFactory();
        final File tempFile;
        final ErrorCheckingPrintWriter tempOut;
        try {
           tempFile = IOUtil.newTempFile("MakeTripletDge", ".tmp", TMP_DIR.toArray(new File[0]));
           tempOut = new ErrorCheckingPrintWriter(tempStreamFactory.wrapTempOutputStream(Files.newOutputStream(tempFile.toPath()), Defaults.BUFFER_SIZE));
        } catch (IOException ex) {
            throw new RuntimeIOException("Error creating temporary file", ex);
        }

        while (!queue.isEmpty()) {
            final TabularDgeStream tabularDgeStream = queue.poll();
            if (!tabularDgeStream.hasNext()) {
                throw new IllegalStateException("Should not be possible");
            }
            final String gene = tabularDgeStream.peekGene();
            if (!gene.equals(currGene)) {
                genes.add(gene);
                currGene = gene;
            }
            final String[] row = tabularDgeStream.next();
            for (int i = 0; i < row.length; ++i) {
                final String value = row[i];
                if (!value.equals(ZERO)) {
                    tempOut.println((genes.size()) + "\t" + (tabularDgeStream.startColumnIndex + i + 1) + "\t" + value);
                    ++nonZeroCount;
                }
            }
            if (tabularDgeStream.hasNext()) {
                queue.add(tabularDgeStream);
            } else {
                CloserUtil.close(tabularDgeStream);
            }
        }
        tempOut.close();
        maybeWriteNamesFiles(genes, outputBarcodes);
        final OutputStream outputStream = IOUtil.openFileForWriting(OUTPUT);
        try {
            final BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(outputStream));
            writer.write(MatrixMarketConstants.makeHeaderLine(ELEMENT_TYPE));
            writer.newLine();
            writer.write(String.format("%d\t%d\t%d", genes.size(), outputBarcodes.size(), nonZeroCount));
            writer.newLine();
            writer.flush();
            final InputStream tempInputStream = tempStreamFactory.wrapTempInputStream(Files.newInputStream(tempFile.toPath()), Defaults.BUFFER_SIZE);
            tempInputStream.transferTo(outputStream);
            outputStream.close();
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
        CloserUtil.close(tempStreamFactory);
        return 0;
    }

    // Yaml keys organized hierarchically
    static class YamlKeys {
        static final String DGES_KEY = "dges";
        static class DatasetsKeys {
            static final String DGE_KEY = "dge";
            static final String PREFIX_KEY = "prefix";
            static final String BARCODE_LIST_KEY = "barcode_list";
            static final String BARCODE_COLUMN_KEY = "barcode_column";
        }
    }

    static class ParsedYamlDge {
        final File dge;
        final String prefix;
        final File barcodeList;
        final String barcodeColumn;

        ParsedYamlDge(File dge, String prefix, File barcodeList, String barcodeColumn) {
            this.dge = dge;
            this.prefix = prefix;
            this.barcodeList = barcodeList;
            this.barcodeColumn = barcodeColumn;
        }
    }

    private Map<String, Object> loadManifest(final File manifest) {
        final Yaml yaml = new Yaml();
        return yaml.load(IOUtil.openFileForReading(manifest));
    }
    private List<ParsedYamlDge> parseYamlDges(Map<String, Object> manifestMap, File manifest) {
        final List<Map<String, Object>> dges = (List<Map<String, Object>>)manifestMap.get(YamlKeys.DGES_KEY);
        if (dges == null) {
            throw new IllegalArgumentException("MANIFEST must contain a '" + YamlKeys.DGES_KEY + "' key");
        }
        File manifestDirectory = manifest.getParentFile();
        if (manifestDirectory == null) {
            manifestDirectory = new File(".");
        }
        final Path manifestDirectoryPath = manifestDirectory.toPath();
        final List<ParsedYamlDge> parsedDges = new ArrayList<>(dges.size());
        for (Map<String, Object> dgeMap : dges) {
            final String dgePath = (String)dgeMap.get(YamlKeys.DatasetsKeys.DGE_KEY);
            if (dgePath == null) {
                throw new IllegalArgumentException("Each dataset must contain a '" + YamlKeys.DatasetsKeys.DGE_KEY + "' key");
            }
            final File dge = manifestDirectoryPath.resolve(dgePath).toFile();
            final String prefix = (String)dgeMap.get(YamlKeys.DatasetsKeys.PREFIX_KEY);
            final String barcodeListPath = (String)dgeMap.get(YamlKeys.DatasetsKeys.BARCODE_LIST_KEY);
            final File barcodeList = barcodeListPath == null ? null : manifestDirectoryPath.resolve(barcodeListPath).toFile();
            final String barcodeColumn = (String)dgeMap.get(YamlKeys.DatasetsKeys.BARCODE_COLUMN_KEY);
            parsedDges.add(new ParsedYamlDge(dge, prefix, barcodeList, barcodeColumn));
        }
        return parsedDges;
    }

    class TabularDgeStreamComparator implements Comparator<TabularDgeStream> {
        @Override
        public int compare(TabularDgeStream o1, TabularDgeStream o2) {
            int ret = o1.peekGene().compareTo(o2.peekGene());
            if (ret == 0) {
                ret = o1.ordinal - o2.ordinal;
            }
            return ret;
        }
    }

    class TabularDgeStream
    implements CloseableIterator<String []> {
        private final File dgeFile;
        private final int ordinal;
        // The number of CBCs from DGEs before this one, so the zero-based index of the first CBC from this DGE
        private final int startColumnIndex;
        private final TabbedTextFileWithHeaderParser parser;
        private final PeekableIterator<TabbedTextFileWithHeaderParser.Row> rowIterator;
        private final List<String> prefixedBarcodes;
        // If all barcodes are included, outputColumnIndices is null.
        private final int[] outputColumnIndices;
        // Recycle the same storage for each row.
        private final String[] filteredValuesForCurrentRow;
        private String prevGene = "";


        /**
         * Prepare to read a DGE file, possibly filtering by a list of barcodes.
         *
         * @param dgeFile the tabular DGE to be read.
         * @param prefix the prefix to prepend to each CBC, along with PREFIX option, or null if no prefix is to be used.
         * @param barcodeList a file containing list of CBCs to include, or null if all CBCs are to be included.
         * @param barcodeColumn if non-null, the name of the column in barcodeList containing the CBCs to include.
         * @param ordinal the zero-based index of this DGE in the manifest.  This is used to keep the order of CBCs
         *                consistent when the same gene appears in multiple DGEs
         * @param startColumnIndex the zero-based index of the first CBC from this DGE in the output.  In essence this is
         *                         the number of CBCs from DGEs before this one.
         */
        TabularDgeStream(final File dgeFile, final String prefix, final File barcodeList, final String barcodeColumn,
                         final int ordinal, final int startColumnIndex) {
            this.dgeFile = dgeFile;
            this.ordinal = ordinal;
            this.startColumnIndex = startColumnIndex;
            IOUtil.assertFileIsReadable(dgeFile);
            Set<String> barcodeSet;
            if (barcodeList != null) {
                IOUtil.assertFileIsReadable(barcodeList);
                if (barcodeColumn == null) {
                    // I can't believe how clumsy this syntax is.
                    barcodeSet = StreamSupport.stream(Spliterators.spliteratorUnknownSize(IOUtil.readLines(barcodeList), Spliterator.ORDERED), false).
                            collect(java.util.stream.Collectors.toSet());
                } else {
                    final TabbedTextFileWithHeaderParser barcodeParser = new TabbedTextFileWithHeaderParser(barcodeList);
                    try {
                        if (!barcodeParser.hasColumn(barcodeColumn)) {
                            throw new IllegalArgumentException("Barcode list " + barcodeList + " does not have column " + barcodeColumn);
                        }
                        barcodeSet = StreamSupport.stream(Spliterators.spliteratorUnknownSize(barcodeParser.iterator(), Spliterator.ORDERED), false).
                                map(row -> row.getField(barcodeColumn)).
                                collect(java.util.stream.Collectors.toSet());
                    } finally {
                        CloserUtil.close(barcodeParser);
                    }
                }
            } else {
                barcodeSet = null;
            }
            // Annoyingly, TabbedTextFileWithHeaderParser doesn't have public access to the column names in order.
            final TabbedInputParser headerParser = new TabbedInputParser(false, dgeFile);
            String[] columnLabels = headerParser.next();
            CloserUtil.close(headerParser);
            if (!columnLabels[0].equals("GENE")) {
                throw new IllegalArgumentException("Header line of " + dgeFile + " does not start with 'GENE'");
            }
            final int numBarcodeColumns = columnLabels.length - 1;
            prefixedBarcodes = new ArrayList<>(columnLabels.length - 1);
            final List<Integer> outputColumnIndicesList = barcodeSet == null? null: new ArrayList<>(columnLabels.length - 1);
            for (int i = 1; i < columnLabels.length; ++i) {
                final String barcode = columnLabels[i];
                if (barcodeSet != null && !barcodeSet.contains(barcode)) {
                    continue;
                }
                final String prefixedBarcode = prefix == null ? barcode : prefix + PREFIX_SEPARATOR + barcode;
                prefixedBarcodes.add(prefixedBarcode);
                if (outputColumnIndicesList != null) {
                    outputColumnIndicesList.add(i);
                }
            }
            if (outputColumnIndicesList != null) {
                outputColumnIndices = outputColumnIndicesList.stream().mapToInt(Integer::intValue).toArray();
                filteredValuesForCurrentRow = new String[outputColumnIndices.length];
            } else {
                outputColumnIndices = null;
                filteredValuesForCurrentRow = new String[numBarcodeColumns];
            }

            parser = new TabbedTextFileWithHeaderParser(dgeFile);
            rowIterator = new PeekableIterator<>(parser.iterator());
        }

        /**
         * Save a little memory by discarding the barcodes once they are written to barcodes.tsv.gz
         */
        void clearBarcodes() {
            prefixedBarcodes.clear();
        }
        public void close() {
            CloserUtil.close(rowIterator);
        }

        @Override
        public boolean hasNext() {
            return rowIterator.hasNext();
        }

        @Override
        public String[] next() {
            final String[] row = rowIterator.next().getFields();
            if (row[0].compareTo(prevGene) <= 0) {
                throw new IllegalArgumentException(String.format("Rows in %s are not sorted by gene at line %d, gene %s",
                        dgeFile,parser.getCurrentLineNumber(), row[0]));
            }
            prevGene = row[0];
            if (outputColumnIndices == null) {
                System.arraycopy(row, 1, filteredValuesForCurrentRow, 0, filteredValuesForCurrentRow.length);
            } else {
                for (int i = 0; i < outputColumnIndices.length; ++i) {
                    filteredValuesForCurrentRow[i] = row[outputColumnIndices[i]];
                }
            }
            return filteredValuesForCurrentRow;
        }

        public String peekGene() {
            return rowIterator.peek().getFields()[0];
        }
    }
}
