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

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Holds a map of smaller->larger barcodes that are unambiguously associated
 * Holds a set of smaller barcodes that are ambiguously associated (and should be removed)
 *
 * @author nemesh
 *
 */
public class BottomUpCollapseResult {
	// maps a smaller barcode that unambiguously is related to a larger barcode
	private Map<String,String> barcodeMap;
	// a set of barcodes that map to > 1 larger barcode
	private Set<String> ambiguousBarcodes;

	public BottomUpCollapseResult() {
		this.barcodeMap=new HashMap<>();
		this.ambiguousBarcodes=new HashSet<>();
	}

	public void addPair (final String smallBC, final String largeBC) {
		if (this.barcodeMap.containsKey(smallBC))
			throw new IllegalArgumentException("This small barcode is already in the map!");
		this.barcodeMap.put(smallBC, largeBC);
	}

	public void addAmbiguousBarcode (final String barcode) {
		this.ambiguousBarcodes.add(barcode);
	}

	public String getLargerRelatedBarcode (final String smallBC) {
		return this.barcodeMap.get(smallBC);
	}

	public Set<String> getAmbiguousBarcodes () {
		return this.ambiguousBarcodes;
	}

	public boolean isAmbiguousBarcode (final String barcode) {
		return this.ambiguousBarcodes.contains(barcode);
	}

	public Set<String> getUnambiguousSmallBarcodes () {
		return this.barcodeMap.keySet();
	}

	public int getCountUnambiguousParentBarcodes () {
		Set<String> s = new HashSet<>(barcodeMap.values());
		return s.size();
	}

	/**
	 * Find if there are any barcodes that are both the intended sequence AND are a smaller target that will be merged.
	 * @return
	 */
	public Set<String> getIntendedAndTargetBarcodes () {
		Set<String> keys=new HashSet<>(barcodeMap.keySet());
		Set<String> values = new HashSet<>(barcodeMap.values());
		keys.retainAll(values);
		return keys;
	}

	public void removeIntendedSequences(final Set<String> barcodes) {
		Set<String> keys = new HashSet<> (this.barcodeMap.keySet());
		for (String smallBarcode: keys) {
			// find if the barcode has an intended sequence.
			String intended = this.barcodeMap.get(smallBarcode);
			// keys can no longer exist if they've been swept out.
			if (intended==null)  continue;
			// if the intended sequence should be removed, then remove the key.
			if (barcodes.contains(intended))
				this.barcodeMap.remove(smallBarcode);
		}
	}

}
