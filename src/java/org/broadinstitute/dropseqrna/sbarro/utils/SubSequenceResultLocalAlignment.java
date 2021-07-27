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

import org.biojava.nbio.alignment.template.AlignedSequence;
import org.biojava.nbio.alignment.template.SequencePair;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.broadinstitute.dropseqrna.utils.editdistance.LevenshteinDistance;
import org.broadinstitute.dropseqrna.utils.editdistance.LevenshteinDistanceResult;

public class SubSequenceResultLocalAlignment implements SubSequenceResultI {

	private final SequencePair<DNASequence, NucleotideCompound> alignedPair;

	// cache results since they are expensive to generate
	private String querySequence = null;
	private String targetSequence = null;
	private String subSequence = null;
	private LevenshteinDistanceResult ldr=null;
	private Integer start;
	private Integer end;

	public SubSequenceResultLocalAlignment (final SequencePair<DNASequence, NucleotideCompound> alignedPair) {
		this.alignedPair=alignedPair;
	}

	@Override
	public LevenshteinDistanceResult getEditDistance() {
		// short circuit if no match found at all.
		if (getSubSequence().length()==0 && this.ldr==null)
			this.ldr = LevenshteinDistance.computeLevenshteinDistanceResult(getQuerySequence(), getTargetSequence());
		// compute if needed.
		if (ldr==null) {
			String s1 = getSubSequence();
			String s2 = getTargetSequence();

			LevenshteinDistanceResult r = LevenshteinDistance.computeLevenshteinDistanceResult(s1, s2);
			this.ldr=r;
		}
		return this.ldr;
	}

	@Override
	public String getQuerySequence() {
		if (this.querySequence==null) {
			AlignedSequence<DNASequence, NucleotideCompound> o = this.alignedPair.getAlignedSequence(1);
			this.querySequence=o.getOriginalSequence().getSequenceAsString();
		}
		return this.querySequence;
	}

	@Override
	public String getSubSequence() {
		if (subSequence==null) {
			// short circuit if there's no bases in common.
			if (alignedPair.getLength()==0) {
				this.subSequence="";
				return this.subSequence;
			}

			String q = getQuerySequence();
			String t = getTargetSequence();
			int start=getStart();
			int end = getEnd();
			if (start==-1 || end==-1)
				this.subSequence=""; // if the start or end is negative, the substring can't be found.
			else {
				// substring is 0 based.
				String subSeq = q.substring(start-1, end);
				this.subSequence=subSeq;
			}
		}
		return this.subSequence;
	}

	@Override
	public int getStart() {
		if (this.start!=null) return this.start;


		if (this.alignedPair.getLength()==0) {
			this.start=-1;
			return -1;
		}

		int s1 = alignedPair.getIndexInTargetAt(1); // The first base in the original target sequence (like an anchor) you aligned to
		int s2 = alignedPair.getIndexInQueryAt(1);  // The first base in the query you aligned to.
		this.start = s2-(s1-1);
		if (this.start <1) this.start=1;

		return this.start;

	}

	@Override
	public int getEnd() {
		if (this.end!=null) return this.end;


		if (this.alignedPair.getLength()==0) {
			this.end=-1;
			return this.end;
		}

		int qLen = getTargetSequence().length();
		int len = alignedPair.getLength();
		int s1 = alignedPair.getIndexInTargetAt(len); // The last base in the original target you aligned to
		int s2 = alignedPair.getIndexInQueryAt(len);  // The last base in the query you aligned to.
		// 	if the aligned length is less than the query length, then this is positive and you need to extend the end coordinate this many more bases
		// if negative, you need to search fewer bases.
		int adjustment = qLen-s1;
		s2=s2+adjustment;
		// can't go past the original length!
		int originalLen = this.alignedPair.getAlignedSequence(1).getOriginalSequence().getLength();
		if (s2>originalLen) s2=originalLen;
		this.end=s2;
		return this.end;
	}

	@Override
	public int getMatchLength() {
		return getSubSequence().length();
	}

	@Override
	public String getTargetSequence() {
		if (this.targetSequence==null) {
			AlignedSequence<DNASequence, NucleotideCompound> o = this.alignedPair.getAlignedSequence(2);
			this.targetSequence= o.getOriginalSequence().getSequenceAsString();
		}

		return this.targetSequence;
	}

	private String removeDeletions (final String seq) {
		StringBuilder b = new StringBuilder();
		for (char c: seq.toCharArray())
			if (c!=DELETION)
				b.append(c);
		return b.toString();
	}

}
