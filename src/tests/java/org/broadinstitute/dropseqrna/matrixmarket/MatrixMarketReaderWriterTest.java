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

import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class MatrixMarketReaderWriterTest {

    @Test(dataProvider = "MMRWTDataProvider")
    public void intTest(final int numRows, final int numCols, final List<String> rowNames, final List<String> colNames,
                        final String rowNamesLabel, final String colNamesLabel) throws IOException {
        final int zeroFrequency = 4;
        final int[][] mat = makeIntMatrix(numRows, numCols, -12345, 3, zeroFrequency);
        final File mmFile = File.createTempFile("MatrixMarketReaderWriterTest.", ".txt.gz");
        mmFile.deleteOnExit();
        final int numNonZeroElements = numRows * numCols - numRows * numCols/zeroFrequency;
        final MatrixMarketWriter writer = new MatrixMarketWriter(mmFile, MatrixMarketConstants.ElementType.integer,
                numRows, numCols, numNonZeroElements, rowNames, colNames, rowNamesLabel, colNamesLabel);
        for (int i = 0; i < numRows; ++i) {
            for (int j = 0; j < numCols; ++j) {
                if (mat[i][j] != 0) {
                    writer.writeTriplet(i, j, mat[i][j]);
                }
            }
        }
        writer.close();

        Assert.assertTrue(MatrixMarketReader.isMatrixMarket(mmFile));
        Assert.assertTrue(MatrixMarketReader.isMatrixMarketInteger(mmFile));

        final MatrixMarketReader reader = new MatrixMarketReader(mmFile, rowNamesLabel, colNamesLabel);
        Assert.assertEquals(reader.getNumRows(), numRows);
        Assert.assertEquals(reader.getNumCols(), numCols);
        Assert.assertEquals(reader.getRowNames(), rowNames == null? Collections.emptyList(): rowNames);
        Assert.assertEquals(reader.getColNames(), colNames == null? Collections.emptyList(): colNames);
        for (final MatrixMarketReader.Element element : reader) {
            Assert.assertNotEquals(element.intValue(), 0);
            Assert.assertEquals(element.intValue(), mat[element.row][element.col]);
            mat[element.row][element.col] = 0;
        }
        for (int i = 0; i < numRows; ++i) {
            for (int j = 0; j < numCols; ++j) {
                if (mat[i][j] != 0) {
                    Assert.fail(String.format("Did not retrieve value at %d, %d", i, j));
                }
            }
        }
        reader.close();
    }

    @Test(dataProvider = "MMRWTDataProvider")
    public void realTest(final int numRows, final int numCols, final List<String> rowNames, final List<String> colNames,
                        final String rowNamesLabel, final String colNamesLabel) throws IOException {
        final int zeroFrequency = 4;
        final double[][] mat = makeDoubleMatrix(numRows, numCols, -12345.6, 3.1, zeroFrequency);
        final File mmFile = File.createTempFile("MatrixMarketReaderWriterTest.", ".txt.gz");
        mmFile.deleteOnExit();
        final int numNonZeroElements = numRows * numCols - numRows * numCols/zeroFrequency;
        final MatrixMarketWriter writer = new MatrixMarketWriter(mmFile, MatrixMarketConstants.ElementType.real,
                numRows, numCols, numNonZeroElements, rowNames, colNames, rowNamesLabel, colNamesLabel);
        for (int i = 0; i < numRows; ++i) {
            for (int j = 0; j < numCols; ++j) {
                if (mat[i][j] != 0) {
                    writer.writeTriplet(i, j, mat[i][j]);
                }
            }
        }
        writer.close();

        Assert.assertTrue(MatrixMarketReader.isMatrixMarket(mmFile));
        Assert.assertTrue(MatrixMarketReader.isMatrixMarketReal(mmFile));

        final MatrixMarketReader reader = new MatrixMarketReader(mmFile, rowNamesLabel, colNamesLabel);
        Assert.assertEquals(reader.getNumRows(), numRows);
        Assert.assertEquals(reader.getNumCols(), numCols);
        Assert.assertEquals(reader.getRowNames(), rowNames == null? Collections.emptyList(): rowNames);
        Assert.assertEquals(reader.getColNames(), colNames == null? Collections.emptyList(): colNames);
        for (final MatrixMarketReader.Element element : reader) {
            Assert.assertNotEquals(element.realValue(), 0);
            Assert.assertNotEquals(mat[element.row][element.col], 0);
            // Check that the ratio of these is close to 1
            Assert.assertEquals(Math.abs(element.realValue()/mat[element.row][element.col]),1.0, 0.000001);
            mat[element.row][element.col] = 0;
        }
        for (int i = 0; i < numRows; ++i) {
            for (int j = 0; j < numCols; ++j) {
                if (mat[i][j] != 0) {
                    Assert.fail(String.format("Did not retrieve value at %d, %d", i, j));
                }
            }
        }
        reader.close();
    }

    @DataProvider(name = "MMRWTDataProvider")
    public Object[][] MMRWTDataProvider() {
        return new Object[][] {
                // Test boundary conditions for number of elements in header lists
                {1000, 1999, makeNames("R", 1000), makeNames("C", 1999), MatrixMarketConstants.GENES, MatrixMarketConstants.CELL_BARCODES},
                // Test default names for row and column name lists
                {100, 1999, makeNames("R", 100), makeNames("C", 1999), null, null},
                // Test no row and column names
                {100, 1999, null, null, null, null},
        };
    }

    private List<String> makeNames(final String prefix, final int numElements) {
        final ArrayList<String> ret = new ArrayList<>(numElements);
        for (int i = 0; i < numElements; ++i) {
            ret.add(prefix + i);
        }
        return ret;
    }

    private int[][] makeIntMatrix(final int numRows, final int numCols, final int initialValue, final int increment, final int zeroFrequency) {
        final int[][] ret = new int[numRows][numCols];
        int val = initialValue;
        int numValues = 0;
        for (int row = 0; row < numRows; row++) {
            for (int col = 0; col < numCols; col++) {
                if (numValues++ % zeroFrequency != 0) {
                    ret[row][col] = val;
                    val += increment;
                }
            }
        }
        return ret;
    }

    private double[][] makeDoubleMatrix(final int numRows, final int numCols, final double initialValue, final double increment, final int zeroFrequency) {
        final double[][] ret = new double[numRows][numCols];
        double val = initialValue;
        int numValues = 0;
        for (int row = 0; row < numRows; row++) {
            for (int col = 0; col < numCols; col++) {
                if (numValues++ % zeroFrequency != 0) {
                    ret[row][col] = val;
                    val += initialValue;
                }
            }
        }
        return ret;
    }
}
