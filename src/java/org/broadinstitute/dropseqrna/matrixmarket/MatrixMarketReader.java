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

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.RuntimeIOException;

import java.io.*;
import java.util.*;
import java.util.function.Function;
import java.util.regex.Pattern;

/**
 * Low-level reader for Matrix Market format.
 * Reads integer or real formats.
 * Optionally supports Drop-seq extension which loads row and column names from comments in the header of the file.
 * Note that the elements returned by the iterator will produce both real and integer values.  If the file is in real
 * format and Element.intValue() is called, the real value will be rounded with Math.round.
 * org.broadinstitute.dropseqrna.priv.matrixmarket.MatrixMarketReader#getElementType() may be used to determine the
 * actual format of the file.
 */
public class MatrixMarketReader
        implements Iterable<MatrixMarketReader.Element>, Closeable {

    /**
     * @param reader Position in file is unchanged by this method.
     * @return true if file looks like Matrix Market
     */
    public static boolean isMatrixMarket(final BufferedReader reader) {
        return (MatrixMarketConstants.MM_HEADER_START.equals(readAndReset(reader, MatrixMarketConstants.MM_HEADER_START.length())));
    }

    /**
     * @param reader Position in file is unchanged by this method.
     * @return true if file looks like Matrix Market integer
     */
    public static boolean isMatrixMarketInteger(final BufferedReader reader) {
        return (MatrixMarketConstants.MM_HEADER_INT.equals(readAndReset(reader, MatrixMarketConstants.MM_HEADER_INT.length())));
    }

    /**
     * @param reader Position in file is unchanged by this method.
     * @return true if file looks like Matrix Market real
     */
    public static boolean isMatrixMarketReal(final BufferedReader reader) {
        return (MatrixMarketConstants.MM_HEADER_REAL.equals(readAndReset(reader, MatrixMarketConstants.MM_HEADER_REAL.length())));
    }

    // Presumably the method references below use the type information to select the appropriate overloaded method
    public static boolean isMatrixMarket(final File file) {
        return isFileTypeHelper(file, MatrixMarketReader::isMatrixMarket);
    }

    public static boolean isMatrixMarketInteger(final File file) {
        return isFileTypeHelper(file, MatrixMarketReader::isMatrixMarketInteger);
    }
    public static boolean isMatrixMarketReal(final File file) {
        return isFileTypeHelper(file, MatrixMarketReader::isMatrixMarketReal);
    }

    private static boolean isFileTypeHelper(final File file, final Function<BufferedReader, Boolean> fun) {
        BufferedReader reader = IOUtil.openFileForBufferedReading(file);
        try {
            return fun.apply(reader);
        } finally {
            CloserUtil.close(reader);
        }
    }

    private static String readAndReset(final BufferedReader reader, final int numChars) {
        try {
            reader.mark(numChars);
            final char[] buf = new char[numChars];
            final int charsRead = reader.read(buf);
            reader.reset();
            if (charsRead != numChars) {
                return null;
            }
            return new String(buf);
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    private final Pattern whitespace = Pattern.compile("\\s+");
    private final String filename;
    private final BufferedReader reader;
    private final int numRows;
    private final int numCols;
    private final int numElements;
    private final List<String> rowNames = new ArrayList<>();
    private final List<String> colNames = new ArrayList<>();
    private final MatrixMarketConstants.ElementType elementType;
    private boolean closed = false;

    /**
     *
     * @param inputFile If has extension .gz, will be opened as a gzip file.
     * @param rowNamesLabel Header label for row names.  If null, "ROWS" is expected.
     * @param colNamesLabel Header label for column names.  If null, "COLS" is expected.
     */
    public MatrixMarketReader(final File inputFile, final String rowNamesLabel, final String colNamesLabel) {
        this(IOUtil.openFileForBufferedReading(inputFile), inputFile.getAbsolutePath(), rowNamesLabel, colNamesLabel);
    }

    public MatrixMarketReader(final File inputFile) {
        this(inputFile, null, null);
    }

    /**
     *
     * @param reader input to be read.
     * @param filename used for error messages only.
     * @param rowNamesLabel Header label for row names.  If null, "ROWS" is expected.
     * @param colNamesLabel Header label for column names.  If null, "COLS" is expected.
     */
    public MatrixMarketReader(final BufferedReader reader, final String filename, String rowNamesLabel, String colNamesLabel) {
        try {
            this.reader = reader;
            this.filename = filename;

            // Determine type of matrix
            final String firstLine = reader.readLine();
            if (MatrixMarketConstants.MM_HEADER_INT.equals(firstLine)) {
                elementType = MatrixMarketConstants.ElementType.integer;
            } else if (MatrixMarketConstants.MM_HEADER_REAL.equals(firstLine)) {
                elementType = MatrixMarketConstants.ElementType.real;
            } else {
                throw new IllegalArgumentException(filename + " does not appear to be a MatrixMarket file.  First line: " +
                        firstLine.substring(0, 200));
            }

            // Read header lists
            if (rowNamesLabel == null) {
                rowNamesLabel = MatrixMarketConstants.ROWS;
            }
            if (colNamesLabel == null) {
                colNamesLabel = MatrixMarketConstants.COLS;
            }
            String line;
            while ((line = reader.readLine()) != null && line.startsWith(MatrixMarketConstants.MM_COMMENT_LINE_START)) {
                if (line.startsWith(MatrixMarketConstants.MM_STRUCTURED_COMMENT_LINE_START)) {
                    line = line.substring(MatrixMarketConstants.MM_STRUCTURED_COMMENT_LINE_START.length());
                    final String[] fields = line.split(MatrixMarketConstants.MM_HEADER_LIST_SEPARATOR);
                    if (fields.length > 1) {
                        if (fields[0].equals(rowNamesLabel)) {
                            rowNames.addAll(Arrays.asList(fields).subList(1, fields.length));
                        } else if (fields[0].equals(colNamesLabel)) {
                            colNames.addAll(Arrays.asList(fields).subList(1, fields.length));
                        }
                    }
                }
            }
            if (line == null) {
                throw new RuntimeException(filename + " appears to be truncated");
            }
            // Get dimensions
            final String[] fields = whitespace.split(line);
            if (fields.length != 3) {
                throw new IllegalArgumentException(filename + " has a matrix dimension line that does not appear to be in MatrixMarket format: " +
                line.substring(0, 200));
            }
            numRows = Integer.parseInt(fields[0]);
            numCols = Integer.parseInt(fields[1]);
            numElements= Integer.parseInt(fields[2]);

            if (rowNames.size() != 0 && rowNames.size() != numRows) {
                throw new RuntimeException(String.format("rowNames.size(%d) != numRows(%d)", rowNames.size(), numRows));
            }
            if (colNames.size() != 0 && colNames.size() != numCols) {
                throw new RuntimeException(String.format("colNames.size(%d) != numCols(%d)", colNames.size(), numCols));
            }
            if ((long)numElements > numRows * (long)numCols) {
                throw new RuntimeException(String.format("numElements(%d) > numRows(%d) * numCols(%d)", numElements, numRows, numCols));
            }

        } catch (IOException e) {
            throw new RuntimeIOException("Exception reading " + filename, e);
        }
    }

    @Override
    public void close() throws IOException {
        if (!closed) {
            closed = true;
            CloserUtil.close(reader);
        }
    }

    public String getFilename() {
        return filename;
    }

    public int getNumRows() {
        return numRows;
    }

    public int getNumCols() {
        return numCols;
    }

    public int getNumElements() {
        return numElements;
    }

    public MatrixMarketConstants.ElementType getElementType() {
        return elementType;
    }

    public List<String> getRowNames() {
        return Collections.unmodifiableList(rowNames);
    }

    public List<String> getColNames() {
        return Collections.unmodifiableList(colNames);
    }

    public static abstract class Element {
        public final int row;
        public final int col;

        Element(int row, int col) {
            this.row = row;
            this.col = col;
        }

        public abstract int intValue();
        public abstract double realValue();
    }

    public static class RealElement extends Element {
        final double val;


        RealElement(int row, int col, double val) {
            super(row, col);
            this.val = val;
        }

        @Override
        public int intValue() {
            return (int)Math.round(val);
        }

        @Override
        public double realValue() {
            return val;
        }
    }

    public static class IntElement extends Element {
        public final int val;

        IntElement(int row, int col, int val) {
            super(row, col);
            this.val = val;
        }

        @Override
        public int intValue() {
            return val;
        }

        @Override
        public double realValue() {
            return val;
        }
    }

    @SuppressWarnings("unchecked")
    @Override
    public Iterator<Element> iterator() {
        if (elementType == MatrixMarketConstants.ElementType.real) {
            return (Iterator) realIterator();
        } else {
            return (Iterator) intIterator();
        }
    }

    public Iterator<RealElement> realIterator() {
        return new MatrixMarketRealIterator();
    }

    public Iterator<IntElement> intIterator() {
        return new MatrixMarketIntIterator();
    }

    private abstract class MatrixMarketIterator {
        private String nextLine;

        MatrixMarketIterator() {
            advance();
        }

        private void advance() {
            try {
                if (!closed) {
                    this.nextLine = reader.readLine();
                    if (this.nextLine == null) {
                        CloserUtil.close(reader);
                        closed = true;
                    }
                } else {
                    this.nextLine = null;
                }
            } catch (IOException e) {
                throw new RuntimeIOException("Exception reading " + filename, e);
            }
        }

        public boolean hasNext() {
            return nextLine != null;
        }

        public Element next() {
            if (!hasNext()) {
                throw new NoSuchElementException();
            }
            final String thisLine = nextLine;
            advance();
            final String fields[] = whitespace.split(thisLine);
            if (fields.length != 3) {
                throw new RuntimeException(filename + " has a bad data line: " + thisLine.substring(0, 200));
            }
            final int oneBasedRow = Integer.parseInt(fields[0]);
            final int oneBasedCol = Integer.parseInt(fields[1]);
            if (oneBasedRow < 1 || oneBasedRow > numRows || oneBasedCol < 1 || oneBasedCol > numCols) {
                throw new RuntimeException("Element line has index out of range: " + thisLine);
            }
            return makeElement(oneBasedRow - 1, oneBasedCol - 1, fields[2]);
        }

        abstract Element makeElement(final int row, final int col, final String val);
    }

    private class MatrixMarketRealIterator
            extends MatrixMarketIterator
            implements Iterator<RealElement> {

        @Override
        public RealElement next() {
            return (RealElement)super.next();
        }

        @Override
        Element makeElement(int row, int col, String val) {
            return new RealElement(row, col, Double.parseDouble(val));
        }
    }

    private class MatrixMarketIntIterator
            extends MatrixMarketIterator
            implements Iterator<IntElement> {

        @Override
        public IntElement next() {
            return (IntElement)super.next();
        }

        @Override
        Element makeElement(int row, int col, String val) {
            return new IntElement(row, col, Integer.parseInt(val));
        }
    }
}
