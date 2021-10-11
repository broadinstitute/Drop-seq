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
package org.broadinstitute.dropseqrna.matrixmarket;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;

import java.io.BufferedWriter;
import java.io.Closeable;
import java.io.File;
import java.io.IOException;
import java.util.List;

public class MatrixMarketWriter
        implements Closeable {

    private final BufferedWriter writer;
    private final String filename;
    private final MatrixMarketConstants.ElementType elementType;
    private final long numRows;
    private final long numCols;
    private final long numNonZeroElements;
    private long numElementsWritten;

    /**
     *
     * @param outputFile If has extension .gz, output file will be gzipped.
     * @param numNonZeroElements the number of triplets to be written
     * @param rowNames If non-null, a list of row names to write into the header.  Length must == numRows
     * @param colNames if non-null, a list of column names to write into the header.  Length must == numCols
     * @param rowNamesLabel If non-null, a label for the rowNames in the header.  If null, "ROWS" is used.
     * @param colNamesLabel If non-null, a label for the rowNames in the header.  If null, "COLS" is used.
     */
    public MatrixMarketWriter(final File outputFile,
                              final MatrixMarketConstants.ElementType elementType,
                              final int numRows,
                              final int numCols,
                              final long numNonZeroElements,
                              final List<String> rowNames,
                              final List<String> colNames,
                              final String rowNamesLabel,
                              final String colNamesLabel) {
        this(IOUtil.openFileForBufferedWriting(outputFile), outputFile.getAbsolutePath(), elementType,
                numRows, numCols, numNonZeroElements, rowNames, colNames, rowNamesLabel, colNamesLabel);
    }

    /**
     *
     * @param filename Used for error messages only
     * @param numNonZeroElements the number of triplets to be written
     * @param rowNames If non-null, a list of row names to write into the header.  Length must == numRows
     * @param colNames if non-null, a list of column names to write into the header.  Length must == numCols
     * @param rowNamesLabel If non-null, a label for the rowNames in the header.  If null, "ROWS" is used.
     * @param colNamesLabel If non-null, a label for the rowNames in the header.  If null, "COLS" is used.
     */
    public MatrixMarketWriter(final BufferedWriter writer,
                              final String filename,
                              final MatrixMarketConstants.ElementType elementType,
                              final int numRows,
                              final int numCols,
                              final long numNonZeroElements,
                              final List<String> rowNames,
                              final List<String> colNames,
                              String rowNamesLabel,
                              String colNamesLabel) {
        try {
            this.writer = writer;
            this.filename = filename;
            this.elementType = elementType;
            this.numRows = numRows;
            this.numCols = numCols;
            this.numNonZeroElements = numNonZeroElements;
            writer.write(MatrixMarketConstants.makeHeaderLine(elementType));
            writer.newLine();
            if (rowNamesLabel == null) {
                rowNamesLabel = MatrixMarketConstants.ROWS;
            }
            if (colNamesLabel == null) {
                colNamesLabel = MatrixMarketConstants.COLS;
            }
            if (rowNames != null) {
                if (numRows != rowNames.size()) {
                    throw new IllegalArgumentException(String.format("numRows(%d) != rownames.size(%d)", numRows, rowNames.size()));
                }
                WriteDimensionLabels(rowNamesLabel, rowNames);
            }
            if (colNames != null) {
                if (this.numCols != colNames.size()) {
                    throw new IllegalArgumentException(String.format("numCols(%d) != colnames.size(%d)", this.numCols, colNames.size()));
                }
                WriteDimensionLabels(colNamesLabel, colNames);
            }
            writer.write(String.format("%d\t%d\t%d", numRows, numCols, numNonZeroElements));
            writer.newLine();
        } catch (IOException e) {
            throw new RuntimeIOException("Exception writing " + filename, e);
        }
    }

    @Override
    public void close() throws IOException {
        writer.close();
    }

    private void WriteDimensionLabels(final String dimensionName, final List<String> names) {
        try {
            for (int i = 0; i < names.size(); i += MatrixMarketConstants.NUM_HEADER_ELEMENTS_PER_ROW) {
                writer.write(MatrixMarketConstants.MM_STRUCTURED_COMMENT_LINE_START);
                writer.write(dimensionName);
                writer.write(MatrixMarketConstants.MM_HEADER_LIST_SEPARATOR);
                final int elementsToWrite = Math.min(MatrixMarketConstants.NUM_HEADER_ELEMENTS_PER_ROW, names.size() - i);
                writer.write(StringUtil.join(MatrixMarketConstants.MM_HEADER_LIST_SEPARATOR, names.subList(i, i + elementsToWrite)));
                writer.newLine();
            }
        } catch (IOException e) {
            throw new RuntimeIOException("Exception writing " + filename, e);
        }
    }

    /**
     * It is legal to call this overload regardless of whether writing integer or real format.
     * @param row 0-based
     * @param col 0-based
     * @param val value to be written
     */
    public void writeTriplet(final long row, final long col, final int val) {
        try {
            assertGoodIndices(row, col);
            writer.write(String.format("%d\t%d\t%d\n", row+1, col+1, val));
        } catch (IOException e) {
            throw new RuntimeIOException("Exception writing " + filename, e);
        }
    }
    /**
     * It is illegal to call this overload if writing integer format.
     * @param row 0-based
     * @param col 0-based
     * @param val value to be written, with 8 digits of precision.
     */
    public void writeTriplet(final long row, final long col, final double val) {
        try {
            if (elementType != MatrixMarketConstants.ElementType.real) {
                throw new UnsupportedOperationException("Cannot write floating-point value to integer matrix");
            }
            assertGoodIndices(row, col);
            writer.write(String.format("%d\t%d\t%.8g\n", row + 1, col + 1, val));
        } catch (IOException e) {
            throw new RuntimeIOException("Exception writing " + filename, e);
        }
    }

    private void assertGoodIndices(final long row, final long col) {
        if (row >= numRows) {
            throw new IllegalArgumentException(String.format("row(%d) >= numRows(%d)", row, numRows));
        }
        if (col >= numCols) {
            throw new IllegalArgumentException(String.format("col(%d) >= numCols(%d)", col, numCols));
        }
        if (numElementsWritten++ >= numNonZeroElements) {
            throw new IllegalArgumentException("More elements written than requested in ctor:" + numNonZeroElements);
        }
    }
}
