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

package org.broadinstitute.dropseqrna.utils.editdistance;

/**
 * Holds information about a pair of barcodes related by Hamming distance of 1.
 * @author nemesh
 *
 */
public class BarcodeSubstitutionPair {

	private final String intendedBarcode;
	private final String neighborBarcode;
	private final int position;
	private final String intendedBase;
	private final String neighborBase;

	public BarcodeSubstitutionPair(final String intendedBarcode, final String neighborBarcode) {
		this.intendedBarcode=intendedBarcode;
		this.neighborBarcode=neighborBarcode;
		int [] pos = HammingDistance.getHammingDistanceChangePositions(neighborBarcode, intendedBarcode);
		if (pos.length!=1)
			new IllegalArgumentException("Strings don't have edit distance 1!");
		this.position=pos[0];
		this.intendedBase = intendedBarcode.substring(pos[0], pos[0]+1);
		this.neighborBase = neighborBarcode.substring(pos[0], pos[0]+1);

	}

	public String getIntendedBarcode() {
		return intendedBarcode;
	}

	public String getNeighborBarcode() {
		return neighborBarcode;
	}

	public int getPosition() {
		return position;
	}

	public String getIntendedBase() {
		return intendedBase;
	}

	public String getNeighborBase() {
		return neighborBase;
	}
}
