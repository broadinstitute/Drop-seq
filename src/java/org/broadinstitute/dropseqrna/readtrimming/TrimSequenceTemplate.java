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
package org.broadinstitute.dropseqrna.readtrimming;



import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.broadinstitute.dropseqrna.TranscriptomeException;

public class TrimSequenceTemplate {

	private final String sequence;
	private final String reverseComplement;
	private final byte[] bases;
	@SuppressWarnings("unused")
	private final byte[] rcBases;
	private byte[] ignoredBases;
	public static final byte A = 'A', C = 'C', G = 'G', T = 'T';

	public TrimSequenceTemplate(String sequence, String ignoredBases) {
		this.sequence = sequence;
		this.reverseComplement = SequenceUtil.reverseComplement(this.sequence);
		bases = StringUtil.stringToBytes(this.sequence);
		rcBases = StringUtil
				.stringToBytes(this.reverseComplement);
		this.ignoredBases = StringUtil
				.stringToBytes(ignoredBases);
	}

	public TrimSequenceTemplate(String barcode) {
		this.sequence = barcode;
		this.reverseComplement = SequenceUtil.reverseComplement(this.sequence);
		bases = StringUtil.stringToBytes(this.sequence);
		rcBases = StringUtil.stringToBytes(this.reverseComplement);
		this.ignoredBases = StringUtil.stringToBytes("Nn");
	}

	public String getSequence() {
		return this.sequence;
	}

	public String reverseComplement() {
		return this.reverseComplement;
	}

	/**
	 * Test to see if this read matches this barcode. If any base of a barcode
	 * starts with N or n, then ignore that position.
	 * 
	 * @param testString
	 *            The read to look for this barcode in. The barcode should be at
	 *            the start of the read for this method.  The entire barcode is expected for a match.
	 * @return true if this barcode is found in the read.
	 */
	public boolean hasForwardMatch(String testString) {
		byte[] testBases = StringUtil.stringToBytes(testString);
		int numBasesCanMatch = 0;
		int numBasesMatch = 0;
		for (int i = 0; i < bases.length; i++) {
			if (isIgnoreBase(this.bases[i]))
				continue;
			numBasesCanMatch++;
			if (SequenceUtil.basesEqual(testBases[i], bases[i])) {
				numBasesMatch++;
			}
		}
		if (numBasesCanMatch == numBasesMatch)
			return (true);
		return false;
	}
	
	/**
	 * Does the testString have part (or all) of the template in it?
	 * The test string starting at position 1 may have some portion of the template.
	 * For example, the last 10 bases of the template would be the first 10 bases of the testString.
	 * The behavior of this is a bit more restricted than getPositionInRead, as it requires you to start finding a match between the template and read,
	 * then continue it to the end of the read.  This is designed to find a some portion of the template at the end of a read, such as an adapter that's glued
	 * on to the end of a read.
	 * @param testString
	 * @param minMatch
	 * @param mismatchesAllowed
	 * @return The position in the template
	 */
	public int getPositionInTemplate (String testString, int minMatch, int mismatchesAllowed) {
		byte [] read = StringUtil.stringToBytes(testString);
		// If the read's too short we can't possibly match it
        if (read == null || read.length < minMatch) return -1;
        
        // read = 48 length
        // template seq = 41 length		
        // Walk forwards
        READ_LOOP:
        for (int start = 0; start<this.bases.length; start++) {
        	final int length = Math.min((this.bases.length - start), read.length);
            int mismatches = 0;
            
            for (int i = 0; i < length; i++) {
            	int idx=i+start;
            	//char base = (char)this.bases[idx];
            	//char readBase = (char) read[i];
                if (!SequenceUtil.isNoCall(this.bases[i]) && !SequenceUtil.basesEqual(this.bases[idx], read[i])) {
                	if (++mismatches > mismatchesAllowed) continue READ_LOOP;
                	
                }
            }
            // If we got this far without breaking out, then it matches.  Make sure it matches at least the min run of bases to return a result.
            int lengthMatch=this.bases.length-start-mismatches;
            if (lengthMatch<minMatch) return (-1);
            return start;
        }
        return -1;
	}
	
	/**
	 * Does the testString have part (or all) of the template in it?
	 * The test string starting at position 1 may have some portion of the template.
	 * For example, ANY bases of the template could have the 8 bases in the middle of a 30 base testString.
	 * This is more flexible and slow than getPositionInTemplate, as it will let any portion of the template occur anywhere in the read.
	 * @param testString A string to test against the template
	 * @param minMatch The number of bases that must match in both 
	 * @param mismatchesAllowed How many mismatches can there be between the template and the read
	 * @return The position in the read, 0 based.
	 */
	public int getPositionInRead (String testString, int minMatch, int mismatchesAllowed) {
		byte [] read = StringUtil.stringToBytes(testString);
		// If the read's too short we can't possibly match it
        if (read == null || read.length < minMatch) return -1;
        
        int maxNumMatches=0;
        int bestReadStartPos=0;
        
        int lastViableReadIndex=read.length-minMatch+1;
        // Walk forwards
        READ_LOOP:  // walks through the read, looking for the template.
        for (int readIndex = 0; readIndex<lastViableReadIndex; readIndex++) {
            int mismatches = 0;
            int numMatches= 0;
            
            // can only search as far as you have left in the read. 
            final int searchLength = Math.min(this.bases.length, read.length-readIndex);
            if (searchLength<minMatch) break;  // if you don't have enough search space left to match a min number of bases you give up.
            for (int templateIndex = 0; templateIndex < searchLength; templateIndex++) {
            	int tempReadIndex=templateIndex+readIndex;
            	//char templateBase = (char)this.bases[templateIndex];
            	//char readBase = (char) read[tempReadIndex];
                if (SequenceUtil.isNoCall(read[tempReadIndex]) || SequenceUtil.basesEqual(this.bases[templateIndex], read[tempReadIndex])) {
                	if (++mismatches > mismatchesAllowed) continue READ_LOOP;
                } else {
                	numMatches++;
                }
                if (numMatches>maxNumMatches) {
                	maxNumMatches=numMatches;
                	bestReadStartPos=readIndex;
                }
            }
            
        }
        if (maxNumMatches<minMatch) return (-1);
        return bestReadStartPos;        
	}
	
	
	

	/**
	 * Does this read have the barcode anywhere in it? The Barcode would be seen
	 * in the reverse orientation, as it's seen on the other strand. This is an
	 * exact match test, so a barcode should not have an ignored base, but the
	 * actual base as observed on the other strand. This is to filter out
	 * barcodes on short inserts.
	 * 
	 * @param readString
	 *            The read string.
	 * @param barcode
	 *            The barcode in the forward orientation.
	 * @return The index of the 5' most base of the barcode in the read, or -1
	 *         if it doesn't exist. This is 1 based - IE: base 1 is the first
	 *         base.
	 */
	public static int findBarcodeIndexReverseMatch(String readString,String barcode) {
		String bc = SequenceUtil.reverseComplement(barcode);
		int index = readString.indexOf(bc);
		if (index > -1)
			index++;
		return index;
	}

	public boolean isIgnoreBase(Byte base) {
		return (baseInBaseList(base, this.ignoredBases));
	}

	public static boolean baseInBaseList(Byte base, byte[] baseList) {
		for (Byte b : baseList) {
			if (SequenceUtil.basesEqual(b, base)) {
				return true;
			}
		}
		return false;
	}

	public byte[] getIgnoredBases() {
		return this.ignoredBases;
	}

	public static Collection<TrimSequenceTemplate> parseBarcodesFromFile(File f) {
		List<TrimSequenceTemplate> result = new ArrayList<TrimSequenceTemplate>();
		try {
			BufferedReader input = new BufferedReader(new FileReader(f));
			try {
				String line = null; // not declared within while loop
				while ((line = input.readLine()) != null) {
					String[] strLine = line.split("\t");
					TrimSequenceTemplate b = new TrimSequenceTemplate(
							strLine[0]);
					result.add(b);
				}
			} finally {
				input.close();
			}
		} catch (IOException ex) {
			throw new TranscriptomeException("Could not read file: "
					+ f.toString());
		}

		return (result);
	}

	/**
	 * If a barcode has ignore bases, then expand those bases to A/C/G/T.
	 * Otherwise, return the barcode. This is recursive, so multiple ignored
	 * bases will be expanded.
	 * 
	 * @return
	 */
	public static Collection<TrimSequenceTemplate> expandBarcode(
			TrimSequenceTemplate b, byte[] ignoredBases) {
		Collection<TrimSequenceTemplate> result = new ArrayList<TrimSequenceTemplate>();
		result.add(b);
		byte[] bases = StringUtil.stringToBytes(b.getSequence());
		for (int i = 0; i < bases.length; i++) {
			boolean ignoreBaseFound = baseInBaseList(bases[i], ignoredBases);
			if (ignoreBaseFound) {
				result.remove(b);
				bases[i] = A;
				TrimSequenceTemplate newBC = new TrimSequenceTemplate(
						StringUtil.bytesToString(bases),
						StringUtil.bytesToString(ignoredBases));
				Collection<TrimSequenceTemplate> r = expandBarcode(newBC,
						ignoredBases);
				result.addAll(r);

				bases[i] = C;
				newBC = new TrimSequenceTemplate(
						StringUtil.bytesToString(bases),
						StringUtil.bytesToString(ignoredBases));
				r = expandBarcode(newBC, ignoredBases);
				result.addAll(r);

				bases[i] = G;
				newBC = new TrimSequenceTemplate(
						StringUtil.bytesToString(bases),
						StringUtil.bytesToString(ignoredBases));
				r = expandBarcode(newBC, ignoredBases);
				result.addAll(r);

				bases[i] = T;
				newBC = new TrimSequenceTemplate(
						StringUtil.bytesToString(bases),
						StringUtil.bytesToString(ignoredBases));
				r = expandBarcode(newBC, ignoredBases);
				result.addAll(r);
				break; // stop looping

			}

		}
		return (result);
	}

	public String toString() {
		return this.sequence;
	}

}
