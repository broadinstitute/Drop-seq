package org.broadinstitute.dropseqrna.beadsynthesis;

import java.util.Map;

public class BiasedBarcodeCollection {
	private Map<String, BeadSynthesisErrorData> biasedBarcodes;
	int numBarcodesFilteredLowUMIs=0;

	public BiasedBarcodeCollection ( final Map<String, BeadSynthesisErrorData> biasedBarcodes, final int numBarcodesFilteredLowUMIs) {
		this.biasedBarcodes=biasedBarcodes;
		this.numBarcodesFilteredLowUMIs=numBarcodesFilteredLowUMIs;
	}

	public Map<String, BeadSynthesisErrorData> getBiasedBarcodes() {
		return biasedBarcodes;
	}

	public int getNumBarcodesFilteredLowUMIs() {
		return numBarcodesFilteredLowUMIs;
	}

}