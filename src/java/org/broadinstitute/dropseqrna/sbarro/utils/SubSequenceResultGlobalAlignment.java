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

public class SubSequenceResultGlobalAlignment implements SubSequenceResultI {

	private final SequencePair<DNASequence, NucleotideCompound> alignedPair;

	// cache results since they are expensive to generate
	private String querySequence = null;
	private String targetSequence = null;
	private String subSequence = null;
	private LevenshteinDistanceResult ldr=null;


	/**
	 * Store the aligned pair of sequences
	 * This class wraps many convenience methods on the alignment result.
	 * @param alignedPair
	 */
	public SubSequenceResultGlobalAlignment (final SequencePair<DNASequence, NucleotideCompound> alignedPair) {
		this.alignedPair=alignedPair;
	}


	/**
	 * Calculate the edit distance between the substring of the best alignment between the query and target, and the original query sequence.
	 * For example, if you queried sequence "ACGT" in the target sequence "AACCTA", the best match would be "ACCT", and would have an edit distance of 1 to "ACGT".
	 * @return The edit distance result
	 */
	@Override
	public LevenshteinDistanceResult getEditDistance () {
		if (ldr==null) {
			int s = this.getStart();
			int e = this.getEnd();

			String s1 = alignedPair.getAlignedSequence(1).getSubSequence(s, e).getSequenceAsString();
			String s2 = alignedPair.getAlignedSequence(2).getSubSequence(s, e).getSequenceAsString();

			LevenshteinDistanceResult r = LevenshteinDistance.computeLevenshteinDistanceResult(s1, s2);
			this.ldr=r;
		}
		return this.ldr;
	}

	/**
	 * The query sequence is the shorter sub-sequence that you find the location of in the larger target sequence.
	 * This is the original sequence you queried for.
	 * @return
	 */
	@Override
	public String getQuerySequence () {
		if (this.querySequence==null) {
			AlignedSequence<DNASequence, NucleotideCompound> o = this.alignedPair.getAlignedSequence(2);
			this.querySequence=o.getOriginalSequence().getSequenceAsString();
		}
		return this.querySequence;
	}

	/**
	 * The best match in the target to the query sequence.
	 * @return
	 */
	@Override
	public String getSubSequence () {
		if (subSequence==null) {
			int s = alignedPair.getAlignedSequence(2).getStart().getPosition();
			int e = alignedPair.getAlignedSequence(2).getEnd().getPosition();
			AlignedSequence<DNASequence, NucleotideCompound> o = alignedPair.getAlignedSequence(1);
			String subSeq = o.getSubSequence(s, e).getSequenceAsString();
			subSeq=removeDeletions(subSeq);
			this.subSequence=subSeq;
		}
		return this.subSequence;
	}

	private String removeDeletions (final String seq) {
		StringBuilder b = new StringBuilder();
		for (char c: seq.toCharArray())
			if (c!=DELETION)
				b.append(c);
		return b.toString();
	}

	/**
	 * Get the start match position of the subsequence in the target sequence.
	 * @return
	 */
	@Override
	public int getStart() {
		return alignedPair.getAlignedSequence(2).getStart().getPosition();
	}

	/**
	 * Get the end match position of the subsequence in the target sequence.
	 * This has to be sensitive to indels!
	 * @return
	 */
	@Override
	public int getEnd() {
		int start = getStart();
		int length = getSubSequence().length();
		return start+length-1;

	}

	/**
	 * How many bases are in the original target sequence that overlap the query?
	 * If bases are missing in the target sequence compared to the query, this can return fewer bases than expected.
	 * If there are extra bases in the target sequence compared to the query, this can return more bases than expected.
	 * @return
	 */
	@Override
	public int getMatchLength() {
		return this.getSubSequence().length();
	}

	/**
	 * The target sequence that you find a subsection of in a query sequence.
	 * @return
	 */
	@Override
	public String getTargetSequence () {
		if (this.targetSequence==null) {
			AlignedSequence<DNASequence, NucleotideCompound> o = this.alignedPair.getAlignedSequence(1);
			this.targetSequence= o.getOriginalSequence().getSequenceAsString();
		}

		return this.targetSequence;
	}

	@Override
	public String toString() {
		return this.alignedPair.toString()+" ED="+this.getEditDistance().getEditDistance();
	}



}
