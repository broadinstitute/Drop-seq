/*
 * MIT License
 *
 * Copyright 2018 Broad Institute
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

import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.editdistance.LevenshteinDistance;
import org.broadinstitute.dropseqrna.utils.editdistance.LevenshteinDistanceResult;

import java.util.List;
import java.util.Map;

public class IntendedSequenceBuilder {

	private final ObjectCounter<String> umiCounts;
	private final Map<String, Double> umiBias;

	public IntendedSequenceBuilder (final ObjectCounter<String> umiCounts, final Map<String, Double> umiBias) {
		this.umiCounts=umiCounts;
		this.umiBias=umiBias;
	}

	public IntendedSequence build (final String intendedSequence, final BarcodeNeighborGroup neighbors) {
		IntendedSequence result = null;
		String root = neighbors.getRootSequence();

		Integer deletedBasePos=null;
		Character deletedBase=null;
		Integer umiCountsIntended = null;
		Double umiBiasIntended = null;


		// need a non-null intendedSequence to check for deletions.
		if (intendedSequence!=null) {
			umiCountsIntended=umiCounts.getCountForKey(intendedSequence);
			umiBiasIntended = umiBias.get(intendedSequence);

			LevenshteinDistanceResult r= LevenshteinDistance.computeLevenshteinDistanceResult(intendedSequence, root, 1, 1, 2);
			String [] ops  = r.getOperations();
			// any position before the last is D, and last is I.

			// gather up the deleted base and position.
			for (int i=0; i<ops.length-2; i++)
				// if a deletion at some base, or a substitution at the last base with an intended sequence, then you can get deletion base/position/rate.
				if (ops[i].equals("D") && ops[ops.length-1].equals("I")) {
					deletedBasePos=i+1; // the position is one based.
					deletedBase = intendedSequence.charAt(deletedBasePos-1); // accessing the array 0 based.
					break;
				}

			// Special case: substitution at the last base only.
			if (ops[ops.length-1].equals("S")) {
				deletedBasePos=ops.length; // the position is one based.
				deletedBase = intendedSequence.charAt(deletedBasePos-1); // accessing the array 0 based.
			}

		}
		int totalRelatedUMIs =neighbors.getNeighborCellBarcodes().stream().mapToInt(x -> this.umiCounts.getCountForKey(x)).sum();
		double [] neighborUMIs = neighbors.getNeighborCellBarcodes().stream().mapToDouble(x -> this.umiCounts.getCountForKey(x)).toArray();
		double [] neighborUMIBias = neighbors.getNeighborCellBarcodes().stream().mapToDouble(x -> this.umiBias.get(x)).toArray();

		Median m = new Median();
		double medianNeighborUMIBias = m.evaluate(neighborUMIBias);
		double medianRelatedSequenceUMIs = m.evaluate(neighborUMIs);

		List<String> neighborBC = neighbors.getNeighborCellBarcodes();

		result = new IntendedSequence(intendedSequence, neighborBC, deletedBase, deletedBasePos,
				umiCountsIntended, medianRelatedSequenceUMIs, totalRelatedUMIs, umiBiasIntended, medianNeighborUMIBias);

		return result;

	}
}
