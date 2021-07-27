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

public class ExtractBarcodeSequences {

	// the templates for searching for anchors in the reads.
	private final FindSubSequence gfpAnchorSearch;
	private final FindSubSequence cassetteAnchorSeach;
	private final int polyABarcodeLength;

	/**
	 * Extracts rabies and anchor sequences from a read.
	 * @param gfpAnchorSequence The GFP anchor barcode sequence (5' of the first rabies barocde, the stop codon barcode)
	 * @param cassetteAnchorSequence The cassette anchor sequence.
	 * @param polyABarcodeLength The length of the polyA rabies barcode.  This is usually 10 bp.
	 */
	public ExtractBarcodeSequences (final String gfpAnchorSequence, final String cassetteAnchorSequence, final int polyABarcodeLength) {
		this.gfpAnchorSearch = new FindSubSequence(gfpAnchorSequence);
		this.cassetteAnchorSeach = new FindSubSequence(cassetteAnchorSequence);
		this.polyABarcodeLength=polyABarcodeLength;
	}

	public ExtractBarcodeSequences (final String gfpAnchorSequence, final String cassetteAnchorSequence) {
		this(gfpAnchorSequence, cassetteAnchorSequence, 10);
	}

	/**
	 * Find the location of the two anchor sequences in the read, then extract the two parts of the rabies virus barcode
	 * based on that sequence.
	 * @param readSequence
	 * @return
	 */
	public ExtractedSequenceGroup findRabiesBarcode (final String readSequence) {
		//TODO: http://winterbe.com/posts/2015/04/07/java8-concurrency-tutorial-thread-executor-examples/
		//TODO: http://www.javapractices.com/topic/TopicAction.do?Id=247
		//TODO: make these all run in parallel?
		SubSequenceResultI gfpAnchor= gfpAnchorSearch.findSequenceLocalAlignment(readSequence);
		SubSequenceResultI cassetteAnchor= cassetteAnchorSeach.findSequenceLocalAlignment(readSequence);
		ExtractedRabiesBarcode stopCodonBarcode=findStopCodonBarcodePositionally(readSequence, gfpAnchor, cassetteAnchor);
		ExtractedRabiesBarcode polyABarcode = findPolyABarcodePositionally(readSequence, cassetteAnchor, this.polyABarcodeLength);

		ExtractedSequenceGroup g = new ExtractedSequenceGroup(gfpAnchor, cassetteAnchor, stopCodonBarcode, polyABarcode);
		return g;
	}

	ExtractedRabiesBarcode findStopCodonBarcodePositionally (final String readSequence, final SubSequenceResultI gfpAnchor, final SubSequenceResultI cassetteAnchor) {

		int start = gfpAnchor.getEnd()+1;
		int end = cassetteAnchor.getStart()-1;
		// start and end are 1 based, substring is 0 based.
		// if the start > end then the anchors could not be detected reasonably, and this result is an empty string with a -1 coordinate.
		if (start > end) return new ExtractedRabiesBarcode("", -1, -1);
		if (start <1 || end <1) return new ExtractedRabiesBarcode("", -1, -1);
		String barcodeSeq = readSequence.substring(start-1, end);
		ExtractedRabiesBarcode result = new ExtractedRabiesBarcode(barcodeSeq, start, end);
		return result;
	}

	/**
	 * This makes the assumption that the polyA rabies barcode is the 10 bases followng the
	 * @param readSequence
	 * @param gfpAnchor
	 * @param cassetteAnchor
	 * @return
	 */
	ExtractedRabiesBarcode findPolyABarcodePositionally (final String readSequence, final SubSequenceResultI cassetteAnchor, final int barcodeLength) {
		int start = cassetteAnchor.getEnd()+1;
		int end = start+barcodeLength-1;
		if (end > readSequence.length())
			end=readSequence.length();
		// if the start > end then the anchors could not be detected reasonably, and this result is an empty string with a -1 coordinate.
		if (start > end) return new ExtractedRabiesBarcode("", -1, -1);
		if (start <1 || end <1) return new ExtractedRabiesBarcode("", -1, -1);
		String barcodeSeq = readSequence.substring(start-1, end);
		ExtractedRabiesBarcode result = new ExtractedRabiesBarcode(barcodeSeq, start, end);
		return result;
	}


}
