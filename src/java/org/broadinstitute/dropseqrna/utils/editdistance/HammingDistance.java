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

