/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
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
