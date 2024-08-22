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

import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.DigitalExpression;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeader;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeaderCodec;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeaderMerger;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.TabbedInputParser;

import java.io.*;
import java.util.*;

@CommandLineProgramProperties(
        summary = "Merges multiple tabular dense Digital Gene Expression files into a single dense tabular file.  " +
                "It is assumed that the cell barcodes in the inputs are non-overlapping, and that genes are in " +
                "alphabetic order.",
        oneLineSummary = "Merges multiple tabular dense Digital Gene Expression files into a single dense tabular file.",
        programGroup = DropSeq.class
)
public class MergeSplitDges
extends CommandLineProgram {
    private final Log log = Log.getInstance(MergeSplitDges.class);

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The DGE files to merge. These must be " +
            "dense, tabular DGEs with gene rows in alphabetic order.")
    public List<File> INPUT;

    @Argument(doc = "The Prefix to add to every cell barcode so that cell barcodes seen "
            + "in multiple DGEs can be differentiated.  A typical setting would be 'EXP1_', yielding barcodes like EXP1_ACTGACCGTTTG.  " +
            "If this is used, it must be invoked as many times as INPUT.  This is option is typically not necessary " +
            "for split DGEs." )
    public List<String> PREFIX;
    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output dense, tabular DGE file.  " +
            "If filename ends with .gz, will be compressed.")
    public File OUTPUT;

    @Argument(doc = "If a gene is not present for an input DGE, use this value for all cell barcodes in that DGE.")
    public String MISSING_VALUE = "0";

    @Argument(shortName = "H", doc="If true, write a header in the DGE file")
    public boolean OUTPUT_HEADER=true;

    @Argument(doc="How strict to be when merging DGE headers. It is expected that split DGEs will contain the same library, " +
            "thus the default is not to complain.")
    public DgeHeaderMerger.Stringency HEADER_STRINGENCY = DgeHeaderMerger.Stringency.NONE;

    @Override
    protected int doWork() {
        INPUT = FileListParsingUtils.expandFileList(INPUT);
        for (final File f: INPUT)
            IOUtil.assertFileIsReadable(f);
        IOUtil.assertFileIsWritable(OUTPUT);
        if (INPUT.isEmpty()) {
            log.error("No input files specified.");
            return 1;
        }
        if (!PREFIX.isEmpty() && INPUT.size()!=PREFIX.size()) {
            log.error("Must have the same number of INPUT and PREFIX parameters, or else no PREFIXes can be specified");
            return 1;
        }
        BufferedWriter out = new BufferedWriter(new OutputStreamWriter(IOUtil.openFileForWriting(OUTPUT)));
        if (OUTPUT_HEADER) {
            final DgeHeader mergedHeader = DgeHeaderMerger.mergeDgeHeaders(INPUT, Collections.EMPTY_LIST, HEADER_STRINGENCY);
            mergedHeader.addCommand(getCommandLine());
            new DgeHeaderCodec().encode(out, mergedHeader);
        }
        final DgeReader[] readers = openInputs();
        final int numCellBarcodes = Arrays.stream(readers).sequential().map(r -> r.numCellBarcodes).reduce(0, Integer::sum);
        final String[] row = new String[numCellBarcodes + 1];

        // Write header row
        row[0] = DigitalExpression.GENE_COLUMN;
        int i = 1;
        for (final DgeReader reader: readers) {
            System.arraycopy(reader.cellBarcodes, 0, row, i, reader.numCellBarcodes);
            i += reader.numCellBarcodes;
        }
        if (i != numCellBarcodes + 1) {
            throw new IllegalStateException("Unpossible! Expected " + numCellBarcodes + " cell barcodes, but found " + i);
        }
        writeLine(out, row);

        // Write gene rows
        String nextGene;
        while ((nextGene = getNextGene(readers)) != null) {
            row[0] = nextGene;
            i = 1;
            for (final DgeReader reader: readers) {
                if (nextGene.equals(reader.nextGene())) {
                    // current gene is present in this DGE so copy its expression values
                    System.arraycopy(reader.nextRow(), 1, row, i, reader.numCellBarcodes);
                } else {
                    // current gene is not present in this DGE so use missing values
                    System.arraycopy(reader.getMissingValues(), 0,  row, i, reader.numCellBarcodes);
                }
                i += reader.numCellBarcodes;
            }
            writeLine(out, row);
        }
        try {
            out.close();
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
        CloserUtil.close(readers);
        return 0;
    }

    private void writeLine(final BufferedWriter out, final String[] row) {
        try {
            out.write(String.join("\t", row));
            out.newLine();
        } catch (final Exception e) {
            throw new RuntimeIOException("Error writing row: " + String.join("\t", row), e);
        }
    }

    /**
     * Open the inputs files and check that there are no duplicate cell barcodes.
     */
    private DgeReader[] openInputs() {
        final DgeReader[] readers = new DgeReader[INPUT.size()];
        final Set<String> cellBarcodes = new HashSet<>();
        for (int i = 0; i < INPUT.size(); ++i) {
            final String prefix = (PREFIX.isEmpty()? null: PREFIX.get(i));
            readers[i] = new DgeReader(INPUT.get(i), prefix);
            for (final String cellBarcode : readers[i].cellBarcodes) {
                if (cellBarcodes.contains(cellBarcode)) {
                    throw new IllegalArgumentException(cellBarcode + " appears in " + INPUT.get(i) +
                            " and some previous input file.");
                }
                cellBarcodes.add(cellBarcode);
            }
        }
        return readers;
    }

    /**
     * Get the next gene in alphabetic order from the readers.
     * @param readers
     * @return the alphabetically earliest gene from the next gene of each input, or null if all inputs are exhausted.
     */
    private String getNextGene(final DgeReader[] readers) {
        String nextGene = null;
        for (final DgeReader reader: readers) {
            final String thisGene = reader.nextGene();
            if (nextGene == null) {
                nextGene = thisGene;
            } else if (thisGene != null && thisGene.compareTo(nextGene) < 0) {
                nextGene = thisGene;
            }
        }
        return nextGene;
    }

    private class DgeReader
    implements Closeable {
        final File file;
        final TabbedInputParser parser;
        private final PeekableIterator<String[]> iterator;
        final String[] cellBarcodes;
        final int numCellBarcodes;
        private String[] missingValues = null;
        private String _nextGene;

        public DgeReader(final File file, final String prefix) {
            this.file = file;
            parser = new TabbedInputParser(false, file);
            iterator = new PeekableIterator<>(parser.iterator());
            final String[] headerRow = iterator.next();
            if (!DigitalExpression.GENE_COLUMN.equals(headerRow[0])) {
                throw new IllegalArgumentException("Header row of " + file + " should start with '" +
                        DigitalExpression.GENE_COLUMN + "', but instead starts with '" + headerRow[0] + "'.");
            }
            numCellBarcodes = headerRow.length - 1;
            cellBarcodes = Arrays.copyOfRange(headerRow, 1, headerRow.length);
            if (prefix != null) {
                for (int i = 0; i < cellBarcodes.length; ++i) {
                    cellBarcodes[i] = prefix + cellBarcodes[i];
                }
            }
            _nextGene = iterator.hasNext()? iterator.peek()[0]: null;
        }

        String nextGene() {
            return _nextGene;
        }

        String[] nextRow() {
            final String[] ret = iterator.next();
            final String peekGene = iterator.hasNext()? iterator.peek()[0]: null;
            if (peekGene != null && _nextGene != null) {
                if (peekGene.compareTo(_nextGene) <= 0) {
                    throw new IllegalArgumentException(String.format("Genes are not in alphabetic order(%s <= %s) in %s",
                            peekGene, _nextGene, file));
                }
            }
            _nextGene = peekGene;
            return ret;
        }

        String[] getMissingValues() {
            if (missingValues == null) {
                missingValues = new String[numCellBarcodes];
                Arrays.fill(missingValues, MISSING_VALUE);
            }
            return missingValues;
        }

        @Override
        public void close() throws IOException {
            parser.close();
        }
    }
}
