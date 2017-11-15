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
package org.broadinstitute.dropseqrna.utils.editdistance;


public class LevenshteinDistance {

	private static int minimum(int a, int b, int c) {
		return Math.min(Math.min(a, b), c);
	}
	
	public static int getDistance (String str1, String str2) {
		return computeLevenshteinDistanceResult(str1, str2).getEditDistance();
	}
	
	public static int getDistance (String str1,String str2, int deletionCost, int insertionCost, int substitutionCost) {
		return computeLevenshteinDistanceResult(str1, str2, deletionCost, insertionCost, substitutionCost).getEditDistance();
	}
 
	public static LevenshteinDistanceResult computeLevenshteinDistanceResult (String str1,String str2) {
		return (computeLevenshteinDistanceResult(str1, str2, 1, 1, 1));
	}
	
	public static LevenshteinDistanceResult computeLevenshteinDistanceResult (String str1,String str2, int deletionCost, int insertionCost, int substitutionCost) {
		int[][] distance = new int[str1.length() + 1][str2.length() + 1];
 
		for (int i = 0; i <= str1.length(); i++)
			distance[i][0] = i;
		for (int j = 1; j <= str2.length(); j++)
			distance[0][j] = j;
 
		for (int i = 1; i <= str1.length(); i++) {
			for (int j = 1; j <= str2.length(); j++) {
				int iCost=distance[i - 1][j] + deletionCost;
				int dCost=distance[i][j - 1] + insertionCost;
				int sCost=distance[i - 1][j - 1]+ ((str1.charAt(i - 1) == str2.charAt(j - 1)) ? 0 : substitutionCost);
				
				int value = minimum(iCost, dCost, sCost);
				distance[i][j] = value;
			}
		}
		LevenshteinDistanceResult r = new LevenshteinDistanceResult(str1, str2, distance, deletionCost, insertionCost, substitutionCost);
		return (r);    
	}
	
	public static int getIndelSlidingWindowEditDistance(String str1, String str2) {
		return (getIndelSlidingWindowEditDistance(str1, str2, Math.max(str1.length(), str2.length())));
	}
	
	/**
	 * Compute the edit distance between two strings, using Jim's indel window modification.
	 * If the naive edit distance is greater than the threshold, then return the threshold.
	 * Otherwise, do lots of fancy (but expensive) stuff to calculate the edit distance.
	 * 
	 * @param str1
	 * @param str2
	 * @param threshold If naive edit distance is larger than this threshold, then return this threshold instead.
	 * @return
	 */
	public static int getIndelSlidingWindowEditDistance(String str1, String str2, int threshold) {
		int r = computeLevenshteinDistanceResult (str1,str2, 1,1,2).getEditDistanceIndelCorrected();
		return r;
	}
	
	/**
	 * 1) Compute comparison string as above. 
	 * a. Starting from the beginning of the comparison string, add 1 to the edit distance for each "S" up to an "I" or "D"
         *
	 * 2) At the first "I" or "D" base at position N in the string:
	 * a. Increase ED by 1.
	 * b. Clip off N bases from barcode A and N-1 bases from barcode B if "D" and N-1 bases from barcode A and N bases from barcode B if "I";
	 * c. Clip off one base from the end of barcode A if an "I" or one base from the end of barcode B if "D";
	 * d. Re-run adist. Increase ED by 1 for each "S";
	 *	i. Return to 2) and repeat until there are no "I" or "D" bases remaining.
	 * @return 
	 */
	/*
	private static int getIndelSlidingWindowEditDistance (String str1, String str2, int editDistance, int threshold) {
		
		LevenshteinDistanceResult r = computeLevenshteinDistanceResult(str1, str2, 1,1,2);
		//LevenshteinDistanceResult r = computeLevenshteinDistanceResult(str1, str2, 1,1,1);
		
		// because we're modifying the substutition cost to be 2x, the threshold should also be 2x.
		if (r.getEditDistance()>threshold*2) return (threshold);
		String [] operations = r.getOperations ();
		
		int idxD=getFirstIndex(operations, "D");
		int idxI=getFirstIndex(operations, "I");
		int idxMin = Math.min(idxD, idxI);
		
		// if there are no INDELS, count the substitutions and exit.
		if (idxMin==-1) {
			int countSubs = countString(operations, "S");
			editDistance+=countSubs;
			return (editDistance);
		}
		
		//otherwise count the substitutions before the deletion
		String [] operationsTrimmed = Arrays.copyOfRange(operations, 0, idxMin+1);
		int countSubs = countString(operationsTrimmed, "S");
		editDistance+=countSubs;
		
		// trim the inputs to substrings
		String str1Trimmed=null;
		String str2Trimmed=null;
		if (idxD==idxMin) {
			str1Trimmed=str1.substring(idxMin+1, str1.length());
			str2Trimmed=str2.substring(idxMin, str2.length()-1);
		}
		
		if (idxI==idxMin) {
			str1Trimmed=str1.substring(idxMin, str1.length()-1);
			str2Trimmed=str2.substring(idxMin+1, str2.length());
		}
		
		editDistance++;
		
		//continue counting edit distance until only "S" or all matches remain
		return (getIndelSlidingWindowEditDistance(str1Trimmed, str2Trimmed, editDistance, threshold));
		
	
	}
	
	private static int getFirstIndex (String [] data, String test) {
		for (int i=0; i<data.length; i++) {
			if (data[i].equals(test)) return (i);
		}
		return -1;
	}
	
	private static int countString (String [] data,  String test) {
		int count=0;
		for (String s: data) {
			if (s.equals(test)) {
				count++;
			}
		}
		return (count);
	}
	*/
}
	
	
