/*
 * MIT License
 *
 * Copyright 2022 Broad Institute
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

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeader;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeaderCodec;
import org.broadinstitute.dropseqrna.matrixmarket.MatrixMarketConstants;
import org.broadinstitute.dropseqrna.matrixmarket.MatrixMarketReader;
import picard.util.TabbedInputParser;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedList;

/**
 * Holds an integer DGE as loaded from sparse or dense input, in triplet format,
 * before sorting or filtering have been done.  Genes rejected by gene enumerator are excluded.
 */
class RawLoadedDge {

    private static final String GENE = "GENE";

    /** sum of transcripts for a cell */
    int[] rawNumTranscripts;
    /** number of non-zero genes for a cell */
    int[] rawNumGenes;
    /** cell barcodes in the order they appear in the input */
    String[] rawCellBarcode;
    /** list of non-zero expression values. geneIndex is from GeneEnumerator, cellIndex is index into the above arrays.  */
    final Collection<SparseDge.Triplet> rawTriplets = new LinkedList<>();
    DgeHeader header;

    RawLoadedDge(final File input, final GeneEnumerator geneEnumerator) {
        final BufferedInputStream inputStream = new BufferedInputStream(IOUtil.openFileForReading(input));
        if (MatrixMarketReader.isMatrixMarketInteger(input)) {
            loadDropSeqSparseDge(inputStream, input, geneEnumerator);
        } else {
            loadTabularDge(inputStream, input, geneEnumerator);
        }
        CloserUtil.close(inputStream);

    }

    private void loadTabularDge(
            final BufferedInputStream inputStream,
            final File input,
            final GeneEnumerator geneEnumerator) {
        this.header = new DgeHeaderCodec().decode(inputStream, input.getAbsolutePath());

        TabbedInputParser parser = new TabbedInputParser(false, inputStream);
        String[] headers = parser.next();
        if (!headers[0].equals(GENE))
            throw new RuntimeException("Unexpected first word in DGE: '" + headers[0] + "' in file " + input.getAbsolutePath());
        // Initialized to 0 by default
        this.rawNumTranscripts = new int[headers.length - 1];
        this.rawNumGenes = new int[headers.length - 1];
        this.rawCellBarcode = Arrays.copyOfRange(headers, 1, headers.length);

        // Read the file
        int i = 1;
        while (parser.hasNext()) {
            ++i;
            final String[] dgeLine = parser.next();
            if (dgeLine.length != this.rawCellBarcode.length + 1)
                throw new RuntimeException("Unexpected number of cellIndex in file " + input.getAbsolutePath() + "; line " + i);
            final int geneId = geneEnumerator.getGeneIndex(dgeLine[0]);
            if (geneId == -1)
                // E.g. for an MT gene
                continue;
            for (int j = 0; j < dgeLine.length - 1; ++j) {
                final String expressionStr = dgeLine[j+1];
                if (expressionStr.equals("0"))
                    continue;
                final int expression = Integer.parseInt(expressionStr);
                this.rawNumTranscripts[j] += expression;
                ++this.rawNumGenes[j];
                this.rawTriplets.add(new SparseDge.Triplet(geneId, j, expression));
            }
        }
    }

    private void loadDropSeqSparseDge(
            final BufferedInputStream inputStream,
            final File input,
            final GeneEnumerator geneEnumerator) {
        final MatrixMarketReader mmReader = new MatrixMarketReader(new BufferedReader(new InputStreamReader(inputStream)),
                input.getAbsolutePath(), MatrixMarketConstants.GENES, MatrixMarketConstants.CELL_BARCODES);
        this.rawNumTranscripts = new int[mmReader.getNumCols()];
        this.rawNumGenes = new int[mmReader.getNumCols()];
        this.rawCellBarcode = mmReader.getColNames().toArray(new String[0]);
        final String[] genes = mmReader.getRowNames().toArray(new String[0]);
        final int[] geneIndices = new int[genes.length];
        for (int i = 0; i < genes.length; ++i)
            geneIndices[i] = geneEnumerator.getGeneIndex(genes[i]);
        for (final MatrixMarketReader.Element element: mmReader) {
            final int geneId = geneIndices[element.row];
            if (geneId == -1)
                // E.g. for an MT gene
                continue;
            final int expression = ((MatrixMarketReader.IntElement) element).val;
            this.rawNumTranscripts[element.col] += expression;
            ++this.rawNumGenes[element.col];
            this.rawTriplets.add(new SparseDge.Triplet(geneId, element.col, expression));
        }
    }
}
