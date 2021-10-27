package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment;

import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import htsjdk.samtools.util.Log;
import picard.util.TabbedTextFileWithHeaderParser;
import picard.util.TabbedTextFileWithHeaderParser.Row;

public class CellContaminationParser {

	private static final Log log = Log.getInstance(CellContaminationParser.class);
	
	public static Map<String,Double> parseCellContamination (File input, Collection<String> expectedCellBarcodes) {
		double maxError=1;
		TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(input);
		
		// validate that the expected header lines are contained
		Set<String> expected = new HashSet<String>(Arrays.asList("cell_barcode", "frac_contamination"));
		Set<String> h = parser.getColumnNames();
		if (!h.containsAll(expected)) {
			parser.close();
			log.error("Parsing cell contamination map, expected header values" + expected.toString());
			log.error("File header ", h.toString());
			throw new IllegalArgumentException("Cell contamination file has the wrong headers!");
		}
		
		Map<String, Double> result = new HashMap<>();
		
		Iterator<Row> iter = parser.iterator();
		
 		// parse each line into the map
		while (iter.hasNext()) {
			Row r = iter.next();
			String cell=r.getField("cell_barcode");
			Double frac = Double.parseDouble(r.getField("frac_contamination"));
			frac=Double.min(frac, maxError);
			// only add something to the map when we expect / need it to save memory.
			if (expectedCellBarcodes.contains(cell))
				result.put(cell, frac);
		}

		parser.close();
		
		final TreeSet<String> onlyInBest = new TreeSet<>(expectedCellBarcodes);		
		onlyInBest.removeAll(result.keySet());
		if (onlyInBest.size()>0) {
			log.warn("Not all cells have a defined ambient RNA estimate!  These will fall back to the default maximum value if supplied, or no cap if not supplied+"
					+ String.join(", ", onlyInBest));
		}
		
		return result;

	}
}
