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
package org.broadinstitute.dropseqrna.beadsynthesis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator;

import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgramProperties;

@CommandLineProgramProperties(
        usage = "For each cell, gather up all the UMIs.  An error in synthesis will result in the last base of the synthesis being fixed in >80% of the UMIs for that cell, across all genes." +
			"This fixed base is T.  Cluster together sets of cell barcodes with UMI bias that are edit distance 1 apart, where the error occurs at the last base.  Search for cell barcode"
			+ "sequences that at ED=1 (using an indel, not a substitution) from one cell barcode sequence to the new cluster of sequences.  If found, this is the intended sequence that explains"
			+ "this cluster of errors.  If not found, the synthesis error is complete or is a deletion at base 12.  "
			+ "For cell barcodes where this occurs, output the cell barcode in a file, as well as (optionally) pad the cell barcodes with N for the error bases.",
        usageShort = "Detect barcode synthesis errors where the final base of a UMI is fixed across all UMIs of a cell.",
        programGroup = DropSeq.class
)

public class DetectIncompleteBeadSynthesisErrors extends AbstractDetectBeadSynthesisErrors {

	private static final Log log = Log.getInstance(DetectIncompleteBeadSynthesisErrors.class);


	private double UMIBiasThreshold=0.6;

    @Override
	protected int doWork() {
    	UMIIterator iterator = prepareUMIIterator();

		BiasedBarcodeCollection biasedBarcodeCollection = findBiasedBarcodes(iterator);
		Map<String, BeadSynthesisErrorData> errorBarcodesWithPositions = biasedBarcodeCollection.getBiasedBarcodes();
		int numCellsFilteredLowUMIs = biasedBarcodeCollection.getNumBarcodesFilteredLowUMIs();
		// based on the errors with positions, find the biased UMIs and cluster.

		List<BarcodeNeighborGroup> neighborGroups = buildBarcodeNeighborGroups(errorBarcodesWithPositions, this.UMIBiasThreshold);

		return 0;

    }


    private List<BarcodeNeighborGroup> buildBarcodeNeighborGroups (final Map<String, BeadSynthesisErrorData> errorBarcodesWithPositions, final double umiBiasThreshold) {
    	Map<String, BarcodeNeighborGroup> result = new HashMap<>();
    	for (String k: errorBarcodesWithPositions.keySet()) {
    		BeadSynthesisErrorData b = errorBarcodesWithPositions.get(k);
    		int polyTErrorPosition = b.getPolyTErrorPosition(umiBiasThreshold);
    		// there's an error if the base isn't -1.
    		if (polyTErrorPosition!= -1) {
    			int umiLength = b.getBaseLength();
    			int numErrors= umiLength-polyTErrorPosition+1;
    			String cellBCRoot = padCellBarcode(b.getCellBarcode(), polyTErrorPosition, umiLength);
    			BarcodeNeighborGroup bng = result.get(cellBCRoot);
    			// if the neighbor group is null create and add it.
    			if (bng==null)
					bng=new BarcodeNeighborGroup();
    			bng.addNeighbor(b);
    			result.put(cellBCRoot, bng);
    		}

    	}

    	return new ArrayList<> (result.values());
    }

    /** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new DetectIncompleteBeadSynthesisErrors().instanceMain(args));
	}
}
