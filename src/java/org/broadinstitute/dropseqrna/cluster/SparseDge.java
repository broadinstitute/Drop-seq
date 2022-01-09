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

import com.google.common.collect.ImmutableList;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeader;

import java.io.BufferedInputStream;
import java.io.File;
import java.lang.reflect.Array;
import java.util.*;

/**
 * Reads a DGE text file (in tabular or Drop-seq Matrix Market format) and stores it in sparse format.
 * Currently, any DGE header is ignored.
 * Cells are sorted in descending order by size.  Genes are filtered by gene enumerator.
 */
public class SparseDge {
    // If memory becomes an issue (which is likely) consider using FastUtil or Trove
    public static class Triplet {
        /** GeneEnumerator ordinal */
        final int geneIndex;
        int cellIndex;
        final int value;

        Triplet(final int geneIndex, final int cellIndex, final int value) {
            this.geneIndex = geneIndex;
            this.cellIndex = cellIndex;
            this.value = value;
        }
    }

    private final File input;
    private final DgeHeader header;
    private int[] numTranscripts;
    private int[] numGenes;
    private String[] cellBarcode;
    private final Collection<Triplet> triplets;
    private final ArrayList<String> discardedCells = new ArrayList<>();

    /**
     * Load a DGE into a sparse in-memory format, with cells in descending order by size.
     * @param input Either tabular DGE text, or Drop-seq Matrix Market sparse format.  May be gzipped.
     * @param geneEnumerator Genes are assigned indices by this.
     */
    public SparseDge(final File input, final GeneEnumerator geneEnumerator) {
        try {
            this.input = input;
            final BufferedInputStream inputStream = new BufferedInputStream(IOUtil.openFileForReading(input));
            final RawLoadedDge rawLoadedDge = new RawLoadedDge(input, geneEnumerator);
            CloserUtil.close(inputStream);
            header = rawLoadedDge.header;
            triplets = rawLoadedDge.rawTriplets;
            sortByCellSize(rawLoadedDge);
        } catch (Exception e) {
            throw new RuntimeException("Problem reading " + input.getAbsolutePath(), e);
        }
    }

    private void sortByCellSize(final RawLoadedDge rawLoadedDge) {

        // Sort cells by rawNumTranscripts (descending)
        // Need Integer in order to use Arrays.sort() with custom comparator
        final Integer[] indices = new Integer[rawLoadedDge.rawCellBarcode.length];
        for (int i = 0; i < indices.length; ++i)
			indices[i] = i;
        Arrays.sort(indices, (o1, o2) -> Integer.compare(rawLoadedDge.rawNumTranscripts[o2], rawLoadedDge.rawNumTranscripts[o1]));

        // Sort cell-oriented outputs according to new sort order
        this.numTranscripts = new int[rawLoadedDge.rawNumTranscripts.length];
        this.numGenes = new int[rawLoadedDge.rawNumGenes.length];
        this.cellBarcode = new String[rawLoadedDge.rawCellBarcode.length];
        for (int i = 0; i < indices.length; ++i) {
            final int unsortedIndex = indices[i];
            this.numTranscripts[i] = rawLoadedDge.rawNumTranscripts[unsortedIndex];
            this.numGenes[i] = rawLoadedDge.rawNumGenes[unsortedIndex];
            this.cellBarcode[i] = rawLoadedDge.rawCellBarcode[unsortedIndex];
        }

        // Renumber triplet according to new sort order
        final int[] oldToNewCellIndexMapping = new int[indices.length];
        for (int i = 0; i < indices.length; ++i)
			oldToNewCellIndexMapping[indices[i]] = i;

        for (final Triplet t: triplets)
			t.cellIndex = oldToNewCellIndexMapping[t.cellIndex];
    }

    public int getNumCells() {
        return cellBarcode.length;
    }

    public Collection<Triplet> getTriplets() {
        return Collections.unmodifiableCollection(triplets);
    }

    public String getCellBarcode(final int cellIndex) {
        return cellBarcode[cellIndex];
    }

    public List<String> getCellBarcodes() {
        return ImmutableList.copyOf(cellBarcode);
    }

    public void prefixCellBarcodes(final String prefix) {
        for (int i = 0; i < cellBarcode.length; ++i)
			cellBarcode[i] = prefix + cellBarcode[i];
    }

    public DgeHeader getHeader() {
        return header;
    }

    public File getFile() {
        return input;
    }

    public void discardSmallestCells(final int numCellsToKeep) {
        if (numCellsToKeep >= this.cellBarcode.length)
			return;
        numTranscripts = Arrays.copyOfRange(numTranscripts, 0, numCellsToKeep);
        numGenes = Arrays.copyOfRange(numGenes, 0, numCellsToKeep);
        cellBarcode = Arrays.copyOfRange(cellBarcode, 0, numCellsToKeep);

        triplets.removeIf(t -> t.cellIndex >= numCellsToKeep);
    }

    public void retainOnlyTheseCells(final Set<String> cellBarcodesToRetain) {
        final BitSet cellsToDiscard = new BitSet(cellBarcode.length);
        for (int i = 0; i < cellBarcode.length; ++i)
			if (!cellBarcodesToRetain.contains(cellBarcode[i]))
				cellsToDiscard.set(i);
        discardCells(cellsToDiscard, false);
    }

    public void discardCellsWithFewGenes(final int minGenes) {
        final BitSet cellsToDiscard = new BitSet(numGenes.length);
        for (int i = 0; i < numGenes.length; ++i)
			if (numGenes[i] < minGenes)
				cellsToDiscard.set(i);
        discardCells(cellsToDiscard, true);
    }

    public void discardCellsWithFewTranscripts(final int minTranscripts) {
        final BitSet cellsToDiscard = new BitSet(numTranscripts.length);
        for (int i = 0; i < numTranscripts.length; ++i)
			if (numTranscripts[i] < minTranscripts)
				cellsToDiscard.set(i);
        discardCells(cellsToDiscard, true);
    }

    private void discardCells(final BitSet cellsToDiscard, boolean shouldCaptureDiscardedCells) {
        if (!cellsToDiscard.isEmpty()) {
            final int[] cellIndexMap = new int[getNumCells()];
            int newCellIndex = 0;
            for (int i = 0; i < cellIndexMap.length; ++i) {
                if (cellsToDiscard.get(i)) {
                    cellIndexMap[i] = -1;
                } else {
                    cellIndexMap[i] = newCellIndex++;
                }
            }

            // We only capture the cells being discarded due to having insufficient number of genes or transcripts
            if (shouldCaptureDiscardedCells) {
                captureDiscardedCellBarcodes(cellsToDiscard);
            }
            numTranscripts = removeElements(numTranscripts, cellsToDiscard);
            numGenes = removeElements(numGenes, cellsToDiscard);
            cellBarcode = removeElements(cellBarcode, cellsToDiscard);
            final Iterator<Triplet> it = triplets.iterator();
            while (it.hasNext()) {
                final Triplet t = it.next();
                if (cellsToDiscard.get(t.cellIndex))
					it.remove();
				else
					t.cellIndex = cellIndexMap[t.cellIndex];
            }
        }
    }

    private void captureDiscardedCellBarcodes(final BitSet cellsToDiscard) {
        discardedCells.ensureCapacity(discardedCells.size() + cellsToDiscard.cardinality());
        cellsToDiscard.stream().forEach(index -> discardedCells.add(cellBarcode[index]));
    }

    private int[] removeElements(final int[] array, final BitSet cellsToDiscard) {
        final int[] ret = new int[array.length - cellsToDiscard.cardinality()];
        removeElementsHelper(array, ret, cellsToDiscard);
        return ret;
    }

    private String[] removeElements(final String[] array, final BitSet cellsToDiscard) {
        final String[] ret = new String[array.length - cellsToDiscard.cardinality()];
        removeElementsHelper(array, ret, cellsToDiscard);
        return ret;
    }

    private void removeElementsHelper(final Object src, final Object dest, final BitSet cellsToDiscard) {
        int from = 0;
        int to = 0;
        final int srcLength = Array.getLength(src);
        while (true) {
            final int indexToDiscard = cellsToDiscard.nextSetBit(from);
            if (indexToDiscard == -1) {
                // Capture trailing cells to include after last discarded cell
                if (srcLength - from > 0)
                    //noinspection SuspiciousSystemArraycopy
                    System.arraycopy(src, from, dest, to, srcLength - from);
                break;
            }
            final int length = indexToDiscard - from;
            if (length > 0)
				//noinspection SuspiciousSystemArraycopy
                System.arraycopy(src, from, dest, to, length);
            to += length;
            from = indexToDiscard + 1;
        }
    }


    public Collection<String> getDiscardedCells() {
        return Collections.unmodifiableCollection(discardedCells);
    }
}
