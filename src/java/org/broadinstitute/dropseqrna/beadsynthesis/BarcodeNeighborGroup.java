package org.broadinstitute.dropseqrna.beadsynthesis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;



public class BarcodeNeighborGroup {

	private Set<BeadSynthesisErrorData> neighbors;

	public BarcodeNeighborGroup() {
		this.neighbors=new HashSet<>();
	}

	public void addNeighbor (final BeadSynthesisErrorData neighbor) {
		this.neighbors.add(neighbor);
	}

	public List<String> getNeighborCellBarcodes () {
		List<String> result = new ArrayList<>();
		for (BeadSynthesisErrorData bsed: neighbors)
			result.add(bsed.getCellBarcode());
		Collections.sort(result);
		return result;
	}

	@Override
	public String toString () {
		StringBuilder b = new StringBuilder();
		b.append("Neighbors " + getNeighborCellBarcodes().toString() +"");
		return (b.toString());
	}
}
