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
package org.broadinstitute.dropseqrna.utils.editdistance;

public class EditDistanceMappingMetric {

	private final String barcode;
	private final int numMergedBarcodes;
	private final int editDistanceUsed;
	private final int editDistanceDiscovered;
	private final int originalObservations;
	private final int totalObservations;
	private final int [] edList;

	public EditDistanceMappingMetric (final String barcode, final int numMergedBarcodes, final int editDistanceUsed, final int editDistanceDiscovered, final int originalObservations, final int totalObservations, final int [] edList) {
		this.barcode=barcode;
		this.numMergedBarcodes=numMergedBarcodes;
		this.editDistanceDiscovered=editDistanceDiscovered;
		this.editDistanceUsed=editDistanceUsed;
		this.totalObservations=totalObservations;
		this.originalObservations=originalObservations;
		this.edList=edList;
	}

	public String getBarcode() {
		return barcode;
	}

	public int getNumMergedBarcodes() {
		return numMergedBarcodes;
	}

	public int getEditDistanceUsed() {
		return editDistanceUsed;
	}

	public int getEditDistanceDiscovered() {
		return editDistanceDiscovered;
	}

	public int[] getEdList() {
		return edList;
	}

	@Override
	public String toString () {
		return "Barcode [" + this.barcode +"] num merged [" + this.numMergedBarcodes +"] ED used ["+this.editDistanceUsed+"] ED discovered [" + this.editDistanceDiscovered +"]";
	}

	public int getTotalObservations() {
		return totalObservations;
	}

	public int getOriginalObservations() {
		return originalObservations;
	}



}
