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

import java.util.ArrayList;
import java.util.List;

public class HammingDistance {

	public static int getHammingDistance(final String sequence1, final String sequence2) {
	    char[] s1 = sequence1.toCharArray();
	    char[] s2 = sequence2.toCharArray();

	    int shorter = Math.min(s1.length, s2.length);
	    int longest = Math.max(s1.length, s2.length);

	    int result = 0;
	    for (int i=0; i<shorter; i++)
			if (s1[i] != s2[i]) result++;

	    result += longest - shorter;

	    return result;
	}

	/**
	 * The two strings must be of the same length to use this test.
	 * @param sequence1 The first String to compare
	 * @param sequence2 The second String to compare
	 * @return The positions at which the strings differ.  The number of elements in the array is the edit distance.
	 */
	public static int [] getHammingDistanceChangePositions(final String sequence1, final String sequence2) {
		char[] s1 = sequence1.toCharArray();
	    char[] s2 = sequence2.toCharArray();
	    if (s1.length!=s2.length)
	    	throw new IllegalArgumentException("Strings not of equal length!  String one ["+sequence1+"] String two ["+ sequence2+"]");
	    List<Integer> resultList = new ArrayList<>();
	    for (int i=0; i<s1.length; i++)
			if (s1[i] != s2[i]) resultList.add(i);
	    int[] result =  resultList.stream().mapToInt(x -> x).toArray();
	    return (result);
	}

	/*
	private String getRandomString(int length) {
		SecureRandom random = new SecureRandom();
		return new BigInteger(length, random).toString(32);

	}
	*/

	public static boolean greaterThanHammingDistance(final String sequence1, final String sequence2, final int minDistance) {
	    char[] s1 = sequence1.toCharArray();
	    char[] s2 = sequence2.toCharArray();

	    int result = 0;
	    for (int i=0; i<s1.length; i++) {
	        if (s1[i] != s2[i]) result++;
	        if (result>minDistance) return true;
	    }
	    return false;
	}

}

