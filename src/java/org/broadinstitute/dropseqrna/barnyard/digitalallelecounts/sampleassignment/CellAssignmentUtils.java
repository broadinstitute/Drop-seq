package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment;

import java.io.File;
import java.util.Collection;
import java.util.Map;

import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.GdacAlleleFrequencyReader;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

/**
 * A hodge-podge of common reusable code across many donor assignment code.
 * @author nemesh
 *
 */
public class CellAssignmentUtils {

	private static final Log log = Log.getInstance(CellAssignmentUtils.class);
	
	/**
	 * Get the value from a map where the map may be null, or the key missing or value null.
	 * @param map A Hashmap that may be null
	 * @param key The key in the hashmap to fetch a value for
	 * @return If the map is null return null.  Otherwise return the value of the key in the map.
	 */
	public static <K,V> V getNullableValue (Map<K,V> map, K key) {
		if (map==null) return null;
		return map.get(key);
	}
	
	/**
	 * Converts a Double to a string.  If the Double is null, emit "NA".
	 * @param d The Double to convert
	 * @return The conversion result.
	 */
	public static String convertNullToString (Double d) {
		if (d==null) return ("NA");
		return d.toString();
	}
	
	
	/**
	 * Gets the most common alternate allele for a variant.
	 * If there is no alternate allele, return 'N' instead.
	 * @param vc The variant context
	 * @return The alternate allele for this context.
	 */
	public static char getAltAllele (VariantContext vc) {
		Allele alt = vc.getAltAlleleWithHighestAlleleCount();
		char altAllele='N';
		if (alt!=null) {
			byte [] altBases = alt.getBases();
			if (altBases.length>0)
				altAllele=StringUtil.byteToChar(altBases[0]);
		}
		return altAllele; 
	}
	
	/**
	 * Parse the contamination file, validate against the expected cell barcodes for the experiment.
	 * @param contaminationMapFile
	 * @param expectedCellBarcodes
	 * @return A map of cell barcode to contamination rate
	 */	
	public static Map<String, Double> getCellContamination (final File inputFile, final Collection<String> expectedCellBarcodes) {
		if (inputFile==null) return null;
		log.info("Per-cell ambient RNA estimates in use to cap likelihoods");
		// otherwise, parse the file and create the likelihood collection.
		Map<String, Double> map = CellContaminationParser.parseCellContamination(inputFile, expectedCellBarcodes);
		return (map);
	}
	
	// Map<String, Double> map = CellContaminationParser.parseCellContamination(contaminationMapFile, expectedCellBarcodes);
		
	public static Map<Interval, Double> getMinorAlleleFrequencyMap (final File inputFile) {
		if (inputFile==null) return (null);
		log.info("Per SNP allele frequency estimates in use to cap likelhiods");
		GdacAlleleFrequencyReader reader = new GdacAlleleFrequencyReader(inputFile);
		return reader.getUmiAlleleFrequencyMap();
	}
		
	
	

}
