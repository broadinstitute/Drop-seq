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
import htsjdk.samtools.util.RuntimeIOException;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

public class CellSizeWriter {
    private static final String CELL_SIZE_HEADER = "CELL_BARCODE\tNUM_TRANSCRIPTS\n";

    private final File cellSizeFile;
    private final BufferedWriter cellSizeWriter;

    public CellSizeWriter(final File cellSizeFile) {
        try {
            this.cellSizeFile = cellSizeFile;
            cellSizeWriter = IOUtil.openFileForBufferedWriting(cellSizeFile);
            cellSizeWriter.write(CELL_SIZE_HEADER);
        } catch (IOException e) {
            throw new RuntimeIOException("Exception writing " + cellSizeFile.getAbsolutePath(), e);
        }
    }

    public void writeSize(final String cellBarcode, final int size) {
        try {
            cellSizeWriter.write(String.format("%s\t%d\n", cellBarcode, size));
        } catch (IOException e) {
            throw new RuntimeIOException("Exception writing " + cellSizeFile.getAbsolutePath(), e);
        }
    }

    public void close() {
        try {
            cellSizeWriter.close();
        } catch (IOException e) {
            throw new RuntimeIOException("Exception closing " + cellSizeFile.getAbsolutePath(), e);
        }
    }
}
