package org.broadinstitute.dropseqrna.beadsynthesis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.dropseqrna.utils.editdistance.EDUtils;

/**
 * Some of the cell barcodes can strongly resemble the primer.
 * Find all cell barcodes that are within EDIT_DISTANCE<X> of the primer and flag them
 * @author nemesh
 *
 */
public class DetectPrimerInUMI {

	private final String primer;
	
	// cache the substrings so you don't have to calculate each time.
	// if you always alter the length of the comparison, this will not help, but that's not the imagined use pattern.
	private Map<Integer, List<String>> primerSubstrings;
	
	public DetectPrimerInUMI(String primer) {
		this.primer=primer;
		primerSubstrings = new HashMap<Integer, List<String>>();
	}
	
	public boolean isStringInPrimer (String str, int editDistance) {
		List<String> primerSubstrings = getSubstrings(primer, str.length());
		Set<String> matchingPrimerSubstrings = EDUtils.getInstance().getStringsWithinHammingDistance(str, primerSubstrings, editDistance);
		return !matchingPrimerSubstrings.isEmpty();
	}
	
	/**
	 * Get all substrings of the string 
	 * @param length How long is each substring
	 * @return A list of substrings
	 */
	public List<String> getSubstrings (String str, int length) {
		if (primerSubstrings.containsKey(length)) return primerSubstrings.get(length);
		List<String> result = new ArrayList<String>();
		
		for (int startPos=0; startPos<(str.length()-length)+1; startPos++) {
			String r = str.substring(startPos, startPos+length);
			result.add(r);
		}
		this.primerSubstrings.put(length, result);
		return result;
	}
	
	
	
}
