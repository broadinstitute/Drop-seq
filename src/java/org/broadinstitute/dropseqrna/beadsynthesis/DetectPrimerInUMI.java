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
		List<String> primerSubstrings = getSubstrings(str.length());
		Set<String> matchingPrimerSubstrings = EDUtils.getInstance().getStringsWithinHammingDistance(str, primerSubstrings, editDistance);
		return !matchingPrimerSubstrings.isEmpty();
	}
	
	/**
	 * Get all substrings of the string 
	 * @param length How long is each substring
	 * @return A list of substrings
	 */
	public List<String> getSubstrings (int length) {
		if (primerSubstrings.containsKey(length)) return primerSubstrings.get(length);
		List<String> result = new ArrayList<String>();
		
		for (int startPos=0; startPos<(this.primer.length()-length)+1; startPos++) {
			String r = this.primer.substring(startPos, startPos+length);
			result.add(r);
		}
		this.primerSubstrings.put(length, result);
		return result;
	}
	
	
	
}
