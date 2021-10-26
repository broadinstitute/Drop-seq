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

import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.SequenceUtil;
import org.biojava.nbio.core.alignment.template.AlignedSequence;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.broadinstitute.dropseqrna.utils.editdistance.LevenshteinDistance;
import org.broadinstitute.dropseqrna.utils.editdistance.LevenshteinDistanceResult;

public class ConsensusSequence {

	private final double alignmentScore;
	private final SequencePair<DNASequence, NucleotideCompound> localAlignment;
	private final boolean readTwoReverseComplimented;
	private byte[] baseQualityReadOne = null;
	private byte[] baseQualityReadTwo = null;

	private String consensusSequence = null;
	private ConsensusSequenceIndex index;
	private static char MISSING_BASE = '-';

	/**
	 * Set up the initial consensus sequence.
	 *
	 * @param globalAlignment
	 * @param alignmentScore
	 * @param readTwoReverseComplimented
	 */
	ConsensusSequence(final SequencePair<DNASequence, NucleotideCompound> localAlignment, final double alignmentScore,
			final boolean readTwoReverseComplimented) {
		this.alignmentScore = alignmentScore;
		this.localAlignment = localAlignment;
		index = new ConsensusSequenceIndex(this.localAlignment);
		this.readTwoReverseComplimented = readTwoReverseComplimented;
	}

	/**
	 * Add the original base qualities of the reads. If read 2 has been reverse
	 * complimented for the best alignment, this will reverse compliment the
	 * base qualities of read 2. If consensus has been calculated before base
	 * qualites were added, this may change the result.
	 *
	 * @param baseQualityReadOne
	 * @param baseQualityReadTwo
	 */
	public void addReadBaseQualities(final String baseQualityReadOne, final String baseQualityReadTwo) {
		testReadLengthBaseQualityLengthMatch(
				this.localAlignment.getAlignedSequence(1).getOriginalSequence().getSequenceAsString(),
				baseQualityReadOne);
		testReadLengthBaseQualityLengthMatch(
				this.localAlignment.getAlignedSequence(2).getOriginalSequence().getSequenceAsString(),
				baseQualityReadTwo);

		this.baseQualityReadOne = SAMUtils.fastqToPhred(baseQualityReadOne);
		this.baseQualityReadTwo = SAMUtils.fastqToPhred(baseQualityReadTwo);
		if (this.readTwoReverseComplimented)
			SequenceUtil.reverseComplement(this.baseQualityReadTwo);
		consensusSequence = null;
	}

	private void testReadLengthBaseQualityLengthMatch(final String seq, final String baseQuality) {
		if (seq.length() != baseQuality.length())
			throw new IllegalArgumentException(
					"Sequence and base quality string are not the same length:" + seq + " " + baseQuality);
	}

	/**
	 * Get the consensus sequence for this read pair.
	 *
	 * @return
	 */
	public String getConsensusSequence() {
		// if (this.baseQualityReadOne==null && this.baseQualityReadTwo==null)
		// return getConsensusSequenceWithoutBaseQuality();
		return getConsensusSequenceWithBaseQuality();
	}

	public String getOriginalReadOne() {
		return this.localAlignment.getAlignedSequence(1).getOriginalSequence().getSequenceAsString();
	}

	public String getOriginalReadTwo() {
		return this.localAlignment.getAlignedSequence(2).getOriginalSequence().getSequenceAsString();
	}

	/**
	 * Compare bases at a read position. For bases that are missing, accept the
	 * non-missing base. For bases that disagree, accept the higher base
	 * quality.
	 *
	 * @param o1
	 * @param o2
	 * @param baseQuality1
	 * @param baseQuality2
	 * @return
	 */
	char getBestSequence(final char o1, final char o2, final byte baseQuality1, final byte baseQuality2) {
		if (o1 == MISSING_BASE)
			return o2;
		if (o2 == MISSING_BASE)
			return o1;
		if (o1 != o2)
			if (baseQuality1 >= baseQuality2)
				return o1;
			else
				return o2;
		return o2;
	}


	/**
	 * Change to a local alignment strategy.
	 *
	 * @return
	 */
	String getConsensusSequenceWithBaseQuality() {

		String foo = index.toString();
		int consensusLength = index.getConsensusLength();
		char[] result = new char[consensusLength];
		char[] seq1 = this.localAlignment.getOriginalSequences().get(0).getSequenceAsString().toCharArray();
		char[] seq2 = this.localAlignment.getOriginalSequences().get(1).getSequenceAsString().toCharArray();

		for (int i = 1; i < consensusLength + 1; i++) {
			// -1 means the base is missing.
			int idx1 = index.getReadOnePosition(i) - 1;
			int idx2 = index.getReadTwoPosition(i) - 1;
			char s1 = this.MISSING_BASE;
			char s2 = this.MISSING_BASE;
			byte q1 = 0;
			byte q2 = 0;
			if (idx1 != -1 & idx1 < seq1.length) {
				s1 = seq1[idx1];
				q1 = this.baseQualityReadOne[idx1];
			}
			if (idx2 != -1 & idx2 < seq2.length) {
				s2 = seq2[idx2];
				q2 = this.baseQualityReadTwo[idx2];
			}
			char r = getBestSequence(s1, s2, q1, q2);
			result[i - 1] = r;
		}

		return new String(result);
	}

	/**
	 * If the start or end bound of the sequence overlap a gap in the original
	 * read sequence, an empty subsequence is returned.
	 *
	 * @param readNum
	 * @param start
	 * @param end
	 * @return
	 */
	public String getOriginalSequenceAtConsensusLocation(final int readNum, final int start, final int end) {
		int startRead = this.index.getReadPosition(readNum, start);
		int endRead = this.index.getReadPosition(readNum, end);
		// search for next best index if the start or end is 0.
		if (startRead==0)
			startRead=searchIndexForNonZero(readNum, start, true);
		if (endRead==0)
			endRead=searchIndexForNonZero(readNum,end, false);

		if (startRead == 0 || endRead == 0)
			return "";
		AlignedSequence<DNASequence, NucleotideCompound> s = localAlignment.getAlignedSequence(readNum);
		return s.getOriginalSequence().getSubSequence(startRead, endRead).getSequenceAsString();
	}

	/**
	 * If the consensusPosition position yields a base not in the read (mapped
	 * to 0, if the position is in a gap in the read), search left or right for
	 * the first non-zero position.
	 *
	 * @param readNum
	 *            1 for 1st read, 2 for 2nd.
	 * @param consensusPosition
	 * @param searchRight
	 * @return
	 */
	private int searchIndexForNonZero(final int readNum, final int consensusPosition, final boolean searchRight) {
		int result = this.index.getReadPosition(readNum, consensusPosition);
		int offSet = 0;
		while (result == 0) {
			if (searchRight)
				offSet++;
			else
				offSet--;
			int idx = consensusPosition+offSet;
			// if you'd fall off the end of the read, stop.
			if (idx==0 || idx > this.index.getConsensusLength())
				return result;
			result = this.index.getReadPosition(readNum, consensusPosition+offSet);
		}
		return result;
	}

	public LevenshteinDistanceResult getLocalAlignmentEditDistance() {
		String s1 = this.localAlignment.getAlignedSequence(1).getSequenceAsString();
		String s2 = this.localAlignment.getAlignedSequence(2).getSequenceAsString();
		return LevenshteinDistance.computeLevenshteinDistanceResult(s1, s2);
	}

	public int getConsensusSequenceLength() {
		return getConsensusSequence().length();
	}

	public boolean isReadTwoReverseComplimented() {
		return readTwoReverseComplimented;
	}

	@Override
	public String toString() {
		StringBuilder b = new StringBuilder();
		b.append(this.localAlignment.toString());
		b.append("ED=" + this.getLocalAlignmentEditDistance().getEditDistance());
		return b.toString();
	}

}
