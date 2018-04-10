/*
 * MIT License
 *
 * Copyright 2018 Broad Institute
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
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

import org.apache.commons.lang.RandomStringUtils;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

public class GenerateRandomUMIs {

	private double umiBiasThreshold=1;
	private char [] bases = {'A', 'C', 'G', 'T'};
	private Random rand = new Random(1);

	public GenerateRandomUMIs (final double umiBiasThreshold) {
		this.umiBiasThreshold=umiBiasThreshold;
	}

	public Collection<String> getUMICollection(final int numBarcodes, final int numBases, final BeadSynthesisErrorType errorType) {

		switch (errorType) {
			case SYNTH_MISSING_BASE:  return getBiasedUMIs(numBarcodes, numBases);
			case PRIMER: throw new IllegalArgumentException("Unsupported error type");
			case SINGLE_UMI: throw new IllegalArgumentException("Unsupported error type");
			case FIXED_FIRST_BASE: throw new IllegalArgumentException("Unsupported error type");
			case OTHER_ERROR: throw new IllegalArgumentException("Unsupported error type");
			default: return getNoErrorUMIs(numBarcodes, numBases);
		}

	}

	/**
	 * Get a collection of random UMIs that are biased at the last base to T.
	 * @param numBarcodes
	 * @param numBases
	 * @return
	 */
	private List<String> getBiasedUMIs (final int numBarcodes, final int numBases) {
		List<char []> s = getNoErrorCharArray(numBarcodes, numBases);
		// set the last base to be T
		for (char [] a: s)
			a[numBases-1]='T';
		return s.stream().map(x -> String.valueOf(x)).collect(Collectors.toList());
	}

	private List<String> getNoErrorUMIs (final int numBarcodes, final int numBases) {
		List<char []> s = getNoErrorCharArray(numBarcodes, numBases);
		return s.stream().map(x -> String.valueOf(x)).collect(Collectors.toList());
	}

	private List<char []> getNoErrorCharArray (final int numBarcodes, final int numBases) {
		List<char []> result = new ArrayList<>();
		// initialize.
		for (int i=0; i<numBarcodes; i++)
			result.add(new char [numBases]);

		// build the [ith] base of all UMIs
		for (int i=0; i<numBases; i++) {

			List<Character> basesAtPosition = getBalancedBaseDistribution(numBarcodes);
			// copy the bases in the base at [i,j].
			for (int j=0; j<numBarcodes; j++) {
				char [] r = result.get(j);
				r[i]=basesAtPosition.get(j);
			}
		}
		return result;
	}

	private List<Character> getBalancedBaseDistribution (final int numBarcodes) {
		ObjectCounter<Character> result = new ObjectCounter<>();
		int numGenerated=0;
		while (numGenerated<numBarcodes) {
			int index = rand.nextInt(4);
			char b = bases[index];
			result.increment(b);
			boolean biased=testBiasedBase(result);
			if (biased)
				// back out the addition
				result.decrement(b);
			else
				numGenerated++;
		}
		// convert the collection of characters into a randomized collection.
		List<Character> charList = new ArrayList<>();
		for (Character c: result.getKeys()) {
			int count = result.getCountForKey(c);
			charList.addAll(Collections.nCopies(count, c));
		}
		Collections.shuffle(charList);
		return charList;
	}

	/**
	 * Is this collection of base counts biased >= the umiBiasThreshold?
	 * @param result
	 * @return
	 */
	private boolean testBiasedBase (final ObjectCounter<Character> result) {
		for (Character c: result.getKeys()) {
			int count = result.getCountForKey(c);
			int total=result.getTotalCount();
			double freq = (double) count / (double) total;
			if (count>1 && freq>=this.umiBiasThreshold)
				return true;
		}
		return false;

	}

	public String getRandomString (final int numBases) {
		return RandomStringUtils.random(numBases, this.bases);
	}
}
