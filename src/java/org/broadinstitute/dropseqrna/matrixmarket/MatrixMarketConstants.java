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

public class MatrixMarketConstants {
    public static final String MM_COMMENT_LINE_START = "%";
    public static final String MM_STRUCTURED_COMMENT_LINE_START = "%%";
    public static final String MM_HEADER_START = MM_STRUCTURED_COMMENT_LINE_START + "MatrixMarket matrix coordinate";
    public static final String MM_MATRIX_TYPE_GENERAL = "general";  // as opposed to symmetric
    public static final String MM_HEADER_REAL = makeHeaderLine(ElementType.real);
    public static final String MM_HEADER_INT = makeHeaderLine(ElementType.integer);
    // The following definitions are for the Drop-seq flavor of Matrix Market, which may contain row and column
    // names in the header comments.
    public static final String MM_HEADER_LIST_SEPARATOR = "\t";
    public static final String ROWS = "ROWS";
    public static final String COLS = "COLS";
    // Legacy names for header lines with row and column names
    public static final String GENES = "GENES";
    public static final String CELL_BARCODES = "CELL_BARCODES";
    // After the header line, DropSeq Matrix Market files start with one of these
    public static final String DROP_SEQ_MATRIX_MARKET_DETECTOR1 = MM_STRUCTURED_COMMENT_LINE_START + GENES + MM_HEADER_LIST_SEPARATOR;
    public static final String DROP_SEQ_MATRIX_MARKET_DETECTOR2 = MM_STRUCTURED_COMMENT_LINE_START + ROWS + MM_HEADER_LIST_SEPARATOR;
    public static final int NUM_HEADER_ELEMENTS_PER_ROW = 1000;

    /**
     * The two Matrix Market data types
     */
    public enum ElementType{real, integer}

    public static String makeHeaderLine(ElementType elementType) {
        return String.format("%s %s %s", MM_HEADER_START, elementType.name(), MM_MATRIX_TYPE_GENERAL);
    }
}
