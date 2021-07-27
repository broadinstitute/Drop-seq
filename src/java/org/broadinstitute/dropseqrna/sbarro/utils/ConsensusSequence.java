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
import org.biojava.nbio.alignment.template.AlignedSequence;
import org.biojava.nbio.alignment.template.SequencePair;
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
	 * If base quality isn't available, then calculate consensus by trusting the
	 * first half of the first read, and the 2nd half of the second read, since
	 * those are closer to the start of the read where the base qualities are
	 * higher and the called base more reliable.
	 *
	 * @return
	 */

	/*
	 * @Deprecated String getConsensusSequenceWithoutBaseQuality() { if
	 * (this.consensusSequence!=null) return consensusSequence; int size =
	 * globalAlignment.getLength();
	 *
	 * AlignedSequence<DNASequence,NucleotideCompound> alignedSeq1 =
	 * globalAlignment.getAlignedSequence(1);
	 * AlignedSequence<DNASequence,NucleotideCompound> alignedSeq2 =
	 * globalAlignment.getAlignedSequence(2);
	 *
	 * char [] seq1 = alignedSeq1.getSequenceAsString().toCharArray(); char []
	 * seq2 = alignedSeq2.getSequenceAsString().toCharArray();
	 *
	 * assert(seq1.length==seq2.length);
	 *
	 * char [] result = new char [seq1.length]; // get the halfway point of the
	 * first read. int seq1Length =
	 * alignedSeq1.getOriginalSequence().getLength(); int halfWayPoint = (int)
	 * Math.floor( seq1Length / 2d);
	 *
	 * for (int i=1; i<=size; i++) { //Returns the column index within an
	 * alignment corresponding to the given index in the original Sequence. One
	 * based. int oIndex1 = alignedSeq1.getSequenceIndexAt(i); // 0 based
	 * arrays. char r = getBestSequence(seq1[i-1], seq2[i-1], oIndex1,
	 * halfWayPoint); result[i-1]=r; } String cs = new String (result);
	 * this.consensusSequence=cs; return cs; }
	 */

	/**
	 * Pick the best base at this point in the sequence, based on the position
	 * in the read. If the sequence index is less than the halfway point, trust
	 * the first sequence more. If the sequence index is greater than the
	 * halfway point, trust the second sequence more.
	 *
	 * @param o1
	 * @param o2
	 * @param seqIndex
	 * @param halfWayPoint
	 * @return
	 */
	/*
	 * @Deprecated char getBestSequence (final char o1, final char o2, final int
	 * seqIndex, final int halfWayPoint) { if (o1==MISSING_BASE) return o2; if
	 * (o2==MISSING_BASE) return o1; if (seqIndex<halfWayPoint) return o1;
	 * return o2; }
	 */

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
	 * In this case, rely on the base quality to dictate how the consensus is
	 * generated. At each position where there is a base call, if there is no
	 * disagreement or only one read has the position, accept the base. At a
	 * position that disagrees, take the higher base quality.
	 *
	 * @return
	 */
	/*
	 * String getConsensusSequenceWithBaseQuality() { // getSequenceIndexAt
	 * returns an index of 1 for positions where there's no base. // those
	 * positions will be '-' though, so we'll use the other read's info. if
	 * (this.consensusSequence!=null) return consensusSequence; int size =
	 * globalAlignment.getLength();
	 *
	 * AlignedSequence<DNASequence,NucleotideCompound> alignedSeq1 =
	 * globalAlignment.getAlignedSequence(1);
	 * AlignedSequence<DNASequence,NucleotideCompound> alignedSeq2 =
	 * globalAlignment.getAlignedSequence(2);
	 *
	 * char [] seq1 = alignedSeq1.getSequenceAsString().toCharArray(); char []
	 * seq2 = alignedSeq2.getSequenceAsString().toCharArray();
	 *
	 * assert(seq1.length==seq2.length);
	 *
	 * char [] result = new char [seq1.length]; for (int i=1; i<=size; i++) {
	 * //Returns the column index within an alignment corresponding to the given
	 * index in the original Sequence. One based. int oIndex1 =
	 * alignedSeq1.getSequenceIndexAt(i); int oIndex2 =
	 * alignedSeq2.getSequenceIndexAt(i); // regular arrays are 0 based... byte
	 * bq1 = this.baseQualityReadOne[oIndex1-1]; byte bq2 =
	 * this.baseQualityReadTwo[oIndex2-1]; // 0 based arrays. char r =
	 * getBestSequence(seq1[i-1], seq2[i-1], bq1, bq2); result[i-1]=r; } String
	 * cs = new String (result); this.consensusSequence=cs; return cs; }
	 */

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
	/*
	 * // getSequenceIndexAt returns an index of 1 for positions where there's
	 * no base. // those positions will be '-' though, so we'll use the other
	 * read's info. if (this.consensusSequence!=null) return consensusSequence;
	 * int size = localAlignment.getLength();
	 *
	 * AlignedSequence<DNASequence,NucleotideCompound> alignedSeq1 =
	 * localAlignment.getAlignedSequence(1);
	 * AlignedSequence<DNASequence,NucleotideCompound> alignedSeq2 =
	 * localAlignment.getAlignedSequence(2);
	 *
	 * int queryStart = localAlignment.getIndexInQueryAt(1); int targetStart =
	 * localAlignment.getIndexInTargetAt(1); int queryEnd =
	 * localAlignment.getIndexInQueryAt(size); int targetEnd =
	 * localAlignment.getIndexInTargetAt(size); int len1=
	 * alignedSeq1.getOriginalSequence().getLength(); int len2
	 * =alignedSeq2.getOriginalSequence().getLength();
	 *
	 *
	 * // determine which read has the most 5' base. // target has most 5' base.
	 * if (targetStart > queryStart) {
	 *
	 * int basesTargetOnly = targetStart - queryStart; int disagreeBasesStart =
	 * queryStart-1; int alignedBases = this.localAlignment.getLength(); int
	 * disagreeBasesEnd = len2 - targetEnd; int basesQueryOnly = len1 - queryEnd
	 * - disagreeBasesEnd; int totalBases =
	 * basesTargetOnly+disagreeBasesStart+alignedBases+disagreeBasesEnd+
	 * basesQueryOnly; } char [] seq1 =
	 * alignedSeq1.getOriginalSequence().getSequenceAsString().toCharArray();
	 * char [] seq2 =
	 * alignedSeq2.getOriginalSequence().getSequenceAsString().toCharArray();
	 *
	 * assert(seq1.length==seq2.length);
	 *
	 * char [] result = new char [seq1.length]; for (int i=1; i<=size; i++) {
	 * //Returns the column index within an alignment corresponding to the given
	 * index in the original Sequence. One based. int oIndex1 =
	 * alignedSeq1.getSequenceIndexAt(i); int oIndex2 =
	 * alignedSeq2.getSequenceIndexAt(i); // regular arrays are 0 based... byte
	 * bq1 = this.baseQualityReadOne[oIndex1-1]; byte bq2 =
	 * this.baseQualityReadTwo[oIndex2-1]; // 0 based arrays. char r =
	 * getBestSequence(seq1[i-1], seq2[i-1], bq1, bq2); result[i-1]=r; } String
	 * cs = new String (result); this.consensusSequence=cs; return cs; }
	 */

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
