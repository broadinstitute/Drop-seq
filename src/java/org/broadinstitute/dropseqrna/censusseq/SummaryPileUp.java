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
package org.broadinstitute.dropseqrna.censusseq;

public class SummaryPileUp {

	private final String donor;
	private int refCount;
	private int altCount;
	private int numSNPs;
	// if we track number of het and alt SNPs, then the math is:
	// (instead of a flat *2)
	// 2 * (# het sites/(#total sites))
	public SummaryPileUp (final String donor) {
		this.donor=donor;
		this.refCount=0;
		this.altCount=0;
		this.numSNPs=0;
	}

	public void incrementNumSNPs () {
		this.numSNPs++;
	}

	public void incrementRefCount (final int count) {
		refCount+=count;
	}

	public void incrementAltCount (final int count) {
		altCount+=count;
	}

	public String getDonor () {
		return this.donor;
	}

	public int getRefCount () {
		return this.refCount;
	}

	public int getAltCount () {
		return this.altCount;
	}

	public int getNumSNPs() {
		return this.numSNPs;
	}

	public double getRatioByAltAlleleFraction () {
		double alt = this.altCount;
		double ref = this.refCount;
		// this is 2*alt because we're considering hets only.  We should probably consider hets and alt/alt sites....
		double result = 2*alt/(alt+ref);
		return (result);
	}


}
