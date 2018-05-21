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

package org.broadinstitute.dropseqrna.utils.editdistance;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang.builder.CompareToBuilder;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

/**
 * How often is each possible base change at a position observed?
 * IE: how often is there an A->C change at base 2.
 *
 * @author nemesh
 *
 */
public class BarcodeSubstitutionCollection {

	private Map<Integer, ObjectCounter<BarcodeSubstitutionElement>> substitutionCountsByPosition;

	public BarcodeSubstitutionCollection() {
		substitutionCountsByPosition = new HashMap<>();
	}

	public void add (final String intendedBarcode, final String neighborBarcode) {
		add(new BarcodeSubstitutionPair(intendedBarcode, neighborBarcode));
	}

	public List<Integer> getPositions () {
		List<Integer> result = new ArrayList<>(this.substitutionCountsByPosition.keySet());
		Collections.sort(result);
		return result;
	}

	public void add (final BarcodeSubstitutionPair p) {
		int pos = p.getPosition();
		BarcodeSubstitutionElement e = new BarcodeSubstitutionElement(p.getIntendedBase(), p.getNeighborBase());
		ObjectCounter<BarcodeSubstitutionElement> o = substitutionCountsByPosition.get(pos);
		if (o==null) {
			o= new ObjectCounter<>();
			substitutionCountsByPosition.put(pos, o);
		}
		o.increment(e);
	}

	public int getCount (final String intendedBase, final String neighborBase, final int position) {
		BarcodeSubstitutionElement e = new BarcodeSubstitutionElement(intendedBase, neighborBase);
		return getCount (e, position);
	}

	public int getCount (final BarcodeSubstitutionElement e, final int position) {
		ObjectCounter<BarcodeSubstitutionElement> o = substitutionCountsByPosition.get(position);
		if (o==null) return 0;
		return o.getCountForKey(e);
	}

	/**
	 * Get the most commonly occurring edit substitution change at a position
	 * @param position The base position to check
	 * @return The most common substitution event
	 */
	public BarcodeSubstitutionElement getMostCommonSubsitution (final int position) {
		ObjectCounter<BarcodeSubstitutionElement> o = substitutionCountsByPosition.get(position);
		if (o==null) return null;
		return o.getMode();
	}

	/**
	 * Get the frequency of the most commonly occurring edit substitution change.
	 * @param position The base position to check
	 * @param e The base change from intended sequence to neighbor sequence.
	 * @return the frequency of this substitution event at this position.
	 */
	public double getSubsitutionFrequency (final BarcodeSubstitutionElement e, final int position) {
		int countMode = getCount (e, position);
		int sum = getTotalSubsitutionsAtPosition(position);
		// no divide by 0 please.
		if (sum==0) return 0;
		double freq = (double) countMode / (double) sum;
		return freq;
	}

	/**
	 * Get the total number of substitution events occurred at this position
	 * @param position
	 * @return
	 */
	public int getTotalSubsitutionsAtPosition (final int position) {
		ObjectCounter<BarcodeSubstitutionElement> o = substitutionCountsByPosition.get(position);
		if (o==null) return 0;
		return o.getTotalCount();
	}


	public BarcodeSubstitutionCollection filterToCommonSubstitutionPatterns (final double freqThreshold) {
		BarcodeSubstitutionCollection result = new BarcodeSubstitutionCollection();

		// scan for common problems.
		List<Integer> positions = this.getPositions();
		for (int i=0; i<positions.size(); i++) {
			BarcodeSubstitutionElement e = getMostCommonSubsitution(i);
			double obeservedFreq = getSubsitutionFrequency(e, i);
			// if observed freq is greater, keep.
			if (obeservedFreq>freqThreshold) {
				int count = this.getCount(e, i);
				ObjectCounter<BarcodeSubstitutionElement> ob = new ObjectCounter<>();
				ob.incrementByCount(e, count);
				result.substitutionCountsByPosition.put(i, ob);
			}
		}
		return result;
	}

	/**
	 * Do the intended and neighbor barcodes substitution change match a known pattern?
	 * Useful after filtering data to common patterns.
	 * @param intendedBarcode
	 * @param neighborBarcode
	 * @return
	 */
	public boolean containsPattern(final String intendedBarcode, final String neighborBarcode) {
		BarcodeSubstitutionPair p = new BarcodeSubstitutionPair(intendedBarcode, neighborBarcode);
		BarcodeSubstitutionElement e = new BarcodeSubstitutionElement(p.getIntendedBase(), p.getNeighborBase());
		int pos = p.getPosition();
		int count = this.getCount(e, pos);
		return (count>0);
	}



	/**
	 * I felt guilty about encoding the pair of base changes as two strings with a ":" in the middle.
	 * @author nemesh
	 *
	 */
	public class BarcodeSubstitutionElement implements Comparable<BarcodeSubstitutionElement> {
		private String intendedBase;
		private String neighborbase;

		public BarcodeSubstitutionElement (final String intendedBase, final String neighborBase) {
			this.intendedBase=intendedBase;
			this.neighborbase=neighborBase;
		}

		public String getIntendedBase() {
			return intendedBase;
		}

		public String getNeighborbase() {
			return neighborbase;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + ((intendedBase == null) ? 0 : intendedBase.hashCode());
			result = prime * result + ((neighborbase == null) ? 0 : neighborbase.hashCode());
			return result;
		}

		@Override
		public boolean equals(final Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			BarcodeSubstitutionElement other = (BarcodeSubstitutionElement) obj;
			if (intendedBase == null) {
				if (other.intendedBase != null)
					return false;
			} else if (!intendedBase.equals(other.intendedBase))
				return false;
			if (neighborbase == null) {
				if (other.neighborbase != null)
					return false;
			} else if (!neighborbase.equals(other.neighborbase))
				return false;
			return true;
		}

		@Override
		public int compareTo(final BarcodeSubstitutionElement o) {
			return new CompareToBuilder().append(this.intendedBase, o.getIntendedBase()).append(this.neighborbase, o.getNeighborbase()).toComparison();
		}

		@Override
		public String toString () {
			return this.intendedBase+":"+this.neighborbase;
		}

	}


}
