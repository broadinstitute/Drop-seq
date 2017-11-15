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


public class HammingDistance {

	public static int getHammingDistance(String sequence1, String sequence2) {
	    char[] s1 = sequence1.toCharArray();
	    char[] s2 = sequence2.toCharArray();

	    int shorter = Math.min(s1.length, s2.length);
	    int longest = Math.max(s1.length, s2.length);

	    int result = 0;
	    for (int i=0; i<shorter; i++) {
	        if (s1[i] != s2[i]) result++;
	    }

	    result += longest - shorter;

	    return result;
	}
	
	/*
	private String getRandomString(int length) {
		SecureRandom random = new SecureRandom();
		return new BigInteger(length, random).toString(32);
		 
	}
	*/
	
	public static boolean greaterThanHammingDistance(String sequence1, String sequence2, int minDistance) {
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

