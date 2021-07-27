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
package org.broadinstitute.dropseqrna.sbarro.utils;

/**
 * Represents a single rabies barcode sequence.
 * @author nemesh
 *
 */
public class ExtractedRabiesBarcode {

	private String sequence; // sequence from read
	private int start;  // start coordinate in read.
	private int end;    // end coordinate in read

	/**
	 *
	 * @param sequence The sequence of the rabies barcode from the read
	 * @param start The start location in the original sequence of the rabies barcode
	 * @param end The end location in the original sequence of the rabies barcode
	 */
	public ExtractedRabiesBarcode (final String sequence, final int start, final int end) {
		this.sequence=sequence;
		this.start=start;
		this.end=end;
	}

	public boolean isValid () {
		return (this.start>=0 & this.end>0);
	}


	public String getSequence() {
		return sequence;
	}



	public int getStart() {
		return start;
	}



	public int getEnd() {
		return end;
	}



	@Override
	public String toString () {
		return sequence+" start["+start+"] end [" + end +"]";
	}

}
