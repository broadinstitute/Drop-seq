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

public class ExtractedSequenceGroup {


	private final SubSequenceResultI gfpAnchor;
	private final SubSequenceResultI cassetteAnchor;

	private final ExtractedRabiesBarcode stopCodonBarcode;
	private final ExtractedRabiesBarcode polyABarcode;

	ExtractedSequenceGroup(final SubSequenceResultI gfpAnchor, final SubSequenceResultI cassetteAnchor, final ExtractedRabiesBarcode stopCodonBarcode, final ExtractedRabiesBarcode polyABarcode) {
		this.gfpAnchor=gfpAnchor;
		this.cassetteAnchor=cassetteAnchor;
		this.stopCodonBarcode=stopCodonBarcode;
		this.polyABarcode=polyABarcode;
	}

	public SubSequenceResultI getGfpAnchor() {
		return gfpAnchor;
	}

	public SubSequenceResultI getCassetteAnchor() {
		return cassetteAnchor;
	}

	public ExtractedRabiesBarcode getStopCodonBarcode() {
		return stopCodonBarcode;
	}

	public ExtractedRabiesBarcode getPolyABarcode() {
		return polyABarcode;
	}

	@Override
	public String toString () {
		StringBuilder b = new StringBuilder();
		b.append("Sequence [" +gfpAnchor.getQuerySequence() +"]\n");
		b.append("GFP anchor [" + gfpAnchor.getSubSequence() +"] Loc[" + gfpAnchor.getStart()+"-"+gfpAnchor.getEnd()+ "] len ["+ gfpAnchor.getMatchLength()+"] ED ["+ gfpAnchor.getEditDistance().getEditDistance() +"]\n");
		b.append("Cassette anchor [" + cassetteAnchor.getSubSequence() +"] loc[" + cassetteAnchor.getStart()+"-"+cassetteAnchor.getEnd()+ "] len ["+ cassetteAnchor.getMatchLength()+"] ED ["+ cassetteAnchor.getEditDistance().getEditDistance() +"]\n");
		b.append("Stop Codon Barcode " + stopCodonBarcode.toString()+"\n");
		b.append("PolyA Barcode [" + polyABarcode.toString());
		return b.toString();
	}



}
