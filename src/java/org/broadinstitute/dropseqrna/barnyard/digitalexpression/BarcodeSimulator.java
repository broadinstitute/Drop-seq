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
package org.broadinstitute.dropseqrna.barnyard.digitalexpression;

import org.apache.commons.lang.RandomStringUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

public class BarcodeSimulator {

	private int numBases;
	char [] bases = {'A', 'C', 'G', 'T', 'N'};
	
	Random rand = new Random();
	
	private List<Integer> positions;
	public BarcodeSimulator (int numBases) {
		this.numBases=numBases;
		// populate the possible positions once.
		positions = new ArrayList<Integer>(numBases);
		for (int i=0; i<numBases; i++) {
			positions.add(i);
		}
	}

	public BarcodeSimulator (int numBases, long seed) {
		this(numBases);
		this.setSeed(seed);
	}

	public BarcodeSimulator (int numBases, char [] bases) {
		this(numBases);
		this.bases = bases;
	}

	public void setSeed(long seed) {
		rand.setSeed(seed);
	}
	
	public String getRandomBarcode() {
		return RandomStringUtils.random(this.numBases, 0, this.bases.length, true, false, this.bases, this.rand);
	}
	
	public List<String> getRandomBarcodes (int numBarcodes) {
		List<String> result = new ArrayList<String>(numBarcodes);
		for (int i=0; i<numBarcodes; i++) {
			result.add(getRandomBarcode());
		}
		return (result);
	}

	/**
	 * Substitutions are single base changes.  For example, from AAA to ATA, the second base "A" is substituted by "T".
	 * @param barcode
	 * @param numSubs
	 * @return
	 */
	public String addSubstitutions (String barcode, int numSubs) {
		char [] b = barcode.toCharArray();
		
		int [] positions = getPositionsToChange(barcode.length(), numSubs);
		for (int i=0; i<positions.length; i++) {
			int pos = positions[i];
			char newBase =  getRandomChar(b[pos]);
			b[pos]=newBase;
		}
		String result = new String(b);
		return result;
	}
	
	/**
	 * Adds an insertion to a barcode at a random position.
	 * What's interesting about this is that barcodes always stay the same length, because we read <X> number of bases as a barcode.
	 * So if a barcode has an insertion, then the same number of bases are trimmed of the tail end.
	 * That's why the Levenstein distance metric is not the standard implementation 
	 * (as the inserted string would be longer, giving you matches at the positions that are now missing) 
	 * @param barcode
	 * @param insertionLength
	 * @return
	 */
	public BarcodeSimulatorResult addInsertion (String barcode, int insertionLength, int position) {
		// max insert position is the (lengthBC) - (insertLength-1)
		//
		List<Character> r = getStringAsCharList(barcode);
		char [] insert = getRandomChars(null, insertionLength);
		List<Character> insert2= new ArrayList<Character>();
		// insert the new bases
		for (char c: insert) insert2.add(c);
		r.addAll(position, insert2);
		// trim the bases after the length, to preserve length.
		r=r.subList(0, barcode.length());
		
		String result = getCharListAsString(r);
		BarcodeSimulatorResult a = new BarcodeSimulatorResult(barcode, result);
		a.addChange(getCharListAsString(insert2), position);

		return a;
	}
	
	private String getCharListAsString (List<Character> s) {
		StringBuilder b = new StringBuilder();
		for (char c: s) {
			b.append(c);
		}
		return b.toString();
	}
	
	
	/**
	 * Adds an insertion to a barcode at a defined position.
	 * What's interesting about this is that barcodes always stay the same length, because we read <X> number of bases as a barcode.
	 * So if a barcode has an insertion, then the same number of bases are trimmed of the tail end.
	 * That's why the Levenstein distance metric is not the standard implementation 
	 * (as the inserted string would be longer, giving you matches at the positions that are now missing) 
	 * If the barcode generated is the same as the input barcode and insertionLength >0, then this tries again until it produces a different barcode.
	 * In the case where your starting barcode was "AAA", and you inserted an "A", you'd get the same output as input.
	 * @param barcode
	 * @param insertionLength
	 * @return
	 */
	public BarcodeSimulatorResult addInsertion (String barcode, int insertionLength) {
		// max insert position is the (lengthBC) - (insertLength-1)
		int maxPos= barcode.length() - (insertionLength-1);
		int pos = getPositionToChange(maxPos);
		BarcodeSimulatorResult result = addInsertion (barcode, insertionLength, pos);
		// if you're gonna return a string without a modification (because bad luck, like inserting an "A" into "AAA") try again.
		if (insertionLength>0 && result.modifiedBarcode.equals(barcode)) {
			return addInsertion(barcode, insertionLength);
		}
				
		return result;
	}
	
	
	private List<Character> getStringAsCharList(String s) {
		List<Character> r = new ArrayList<Character>(s.length());
		for (char c: s.toCharArray()) {
			r.add(c);
		}
		return (r);
	}
	
	
	
	
	private char [] getRandomChars(Character excludeCharacter, int count) {
		char [] r = new char [count];
		for (int i=0; i<count; i++) {
			r[i]=getRandomChar(excludeCharacter);
		}
		return (r);
	}
	
	
	/**
	 * Generate a random character that isn't this character.
	 * If excludeCharacter is null, then any of the bases is fine, and is more efficient.
	 * @param excludeCharacter
	 * @return
	 */
	private char getRandomChar(Character excludeCharacter) {
		char [] charSet;
		if (excludeCharacter==null) {
			charSet=this.bases;
			return charSet[getPositionToChange(charSet.length)];
		}
		// if current isn't null, then generate a character list and pick one at random.
		List<Character> result = new ArrayList<Character>(this.bases.length);
		for (char c: this.bases) {
			if (c!=excludeCharacter) {
				result.add(c);
			}
		}
		int index = getPositionToChange(result.size());
		return (result.get(index));
		
	}
	
	
	/**
	 * Functions as a sample without replacement.  Generates multiple positions to change that do not overlap.
	 * 
	 * @param numPositions
	 * @param numChanges
	 * @return
	 */
	private int [] getPositionsToChange (int numPositions, int numChanges) {
		int [] result = new int [numChanges];
		Collections.shuffle(this.positions);
		
		for (int i=0; i<numChanges; i++) {
			result[i]=this.positions.get(i);
		}
		return (result);
	}
	
	/**
	 * This is faster for generating a single position to change.
	 * If called twice, can generate the same position two times.
	 * @return
	 */
	private int getPositionToChange(int numBases) {
		 return rand.nextInt((numBases));
	}
	
	public class BarcodeSimulatorResult {
		
		String originalBarcode;
		String modifiedBarcode;
		List<String> changes;
		List<Integer> positions;
		
		public BarcodeSimulatorResult (String originalBarcode, String modifiedBarcode) {
			this.originalBarcode=originalBarcode;
			this.modifiedBarcode=modifiedBarcode;
			changes = new ArrayList<String>();
			positions = new ArrayList<Integer>();
		}
		
		public void addChange (String change, int position) {
			this.changes.add(change);
			this.positions.add(position);
		}
		
		public List<String> getChanges() {
			return this.changes;
		}
		
		public List<Integer> getPositions() {
			return this.positions;
		}
		
		public String toString () {
			StringBuilder b = new StringBuilder();
			b.append("Original: " + originalBarcode +" ");
			b.append("Modified: " + modifiedBarcode + " ");
			
			for (int i=0; i<changes.size(); i++) {
				b.append("[" + positions.get(i) +"] " +changes.get(i) +" ");
			}
			return b.toString();
		}
		
	}
	
}
