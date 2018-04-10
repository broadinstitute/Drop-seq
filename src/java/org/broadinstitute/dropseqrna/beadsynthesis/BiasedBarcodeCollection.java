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

import java.util.HashMap;
import java.util.Map;

import org.broadinstitute.dropseqrna.utils.ObjectCounter;

public class BiasedBarcodeCollection {
	private final Map<String, BeadSynthesisErrorData> biasedBarcodes;
	private final ObjectCounter<String> umiCount;
	private final Map<String, Double> umiBias;

	public BiasedBarcodeCollection ( final Map<String, BeadSynthesisErrorData> biasedBarcodes) {
		this(biasedBarcodes, new ObjectCounter<String>(), new HashMap<String, Double>());
	}

	public BiasedBarcodeCollection ( final Map<String, BeadSynthesisErrorData> biasedBarcodes, final ObjectCounter<String> allCellBarcodes, final Map<String, Double> umiBias) {
		this.biasedBarcodes=biasedBarcodes;
		this.umiCount=allCellBarcodes;
		this.umiBias=umiBias;
	}

	public void add (final BeadSynthesisErrorData data) {
		this.biasedBarcodes.put(data.getCellBarcode(), data);
	}

	public void addUMIBias (final String cellBarcode, final Double bias) {
		umiBias.put(cellBarcode, bias);
	}

	public void incrementUMICount (final String cellBarcode, final int umiCount) {
		this.umiCount.incrementByCount(cellBarcode, umiCount);
	}

	public Map<String, BeadSynthesisErrorData> getBiasedBarcodes() {
		return biasedBarcodes;
	}

	public ObjectCounter<String> getUMICounts() {
		return umiCount;
	}

	public Map<String, Double> getUMIBias () {
		return this.umiBias;
	}



}
