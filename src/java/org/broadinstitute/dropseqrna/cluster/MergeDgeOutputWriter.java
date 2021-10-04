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

import htsjdk.samtools.util.RuntimeIOException;
import org.broadinstitute.dropseqrna.matrixmarket.MatrixMarketConstants;
import org.broadinstitute.dropseqrna.matrixmarket.MatrixMarketWriter;

import java.io.File;
import java.io.IOException;
import java.util.List;

/**
 * Writer for raw and scaled DGE in Drop-seq Matrix Market format.
 */
class MergeDgeOutputWriter {

    private final MatrixMarketWriter rawDgeWriter;
    private final MatrixMarketWriter scaledDgeWriter;

    /**
     * Prepare to write Drop-seq style Matrix Market sparse matrix files
     * At least one of rawDgeFile and scaledDgeFile must be non-null
     * @param rawDgeFile output for raw (integer) output
     * @param scaledDgeFile output for scaled (columns sum to 1) output
     * @param numNonZeroElements total non-zero elements that will be written
     * @param genes row names.
     * @param cellBarcodes column names.
     */
    public MergeDgeOutputWriter(final File rawDgeFile, final File scaledDgeFile, final long numNonZeroElements,
                                final List<String> genes, final List<String> cellBarcodes) {
        if (rawDgeFile == null && scaledDgeFile == null) {
            throw new IllegalArgumentException("Doesn't make sense to construct with both files null");
        }
        if (rawDgeFile != null) {
            rawDgeWriter = new MatrixMarketWriter(rawDgeFile, MatrixMarketConstants.ElementType.integer,
                    genes.size(), cellBarcodes.size(), numNonZeroElements, genes, cellBarcodes,
                    MatrixMarketConstants.GENES, MatrixMarketConstants.CELL_BARCODES);
        } else {
            rawDgeWriter = null;
        }
        if (scaledDgeFile != null) {
            scaledDgeWriter = new MatrixMarketWriter(scaledDgeFile, MatrixMarketConstants.ElementType.real,
                    genes.size(), cellBarcodes.size(), numNonZeroElements, genes, cellBarcodes,
                    MatrixMarketConstants.GENES, MatrixMarketConstants.CELL_BARCODES);
        } else {
            scaledDgeWriter = null;
        }
    }

    public void close() {
        try {
            if (rawDgeWriter != null) {
                rawDgeWriter.close();
            }
            if (scaledDgeWriter != null) {
                scaledDgeWriter.close();
            }
        } catch (IOException e) {
            throw new RuntimeIOException(e);
        }
    }

    /**
     *
     * @param geneIndex zero-based row index
     * @param cellIndex zero-based column index
     * @param raw       integer DGE value
     * @param scaled    scaled by total DGE for the cell.
     */
    public void writeValue(final int geneIndex, final int cellIndex, final int raw, final double scaled) {
        if (rawDgeWriter != null) {
            rawDgeWriter.writeTriplet(geneIndex, cellIndex, raw);
        }
        if (scaledDgeWriter != null) {
            scaledDgeWriter.writeTriplet(geneIndex, cellIndex, scaled);
        }
    }
}
