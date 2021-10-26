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

import org.biojava.nbio.core.alignment.template.AlignedSequence;
import org.biojava.nbio.core.alignment.template.SequencePair;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;

import java.util.Arrays;

/**
 * A helper class that holds a pair of aligned sequences and figures out all the
 * offsets to the original sequences
 *
 * @author nemesh
 *
 */
public class ConsensusSequenceIndex {

	//private final int firstReadIndex;
	// private final int secondReadIndex;

	private final SequencePair<DNASequence, NucleotideCompound> localAlignment;
	//private final int [][] originalSequencePositions;
	private final int [] [] index;

	public ConsensusSequenceIndex(final SequencePair<DNASequence, NucleotideCompound> localAlignment) {
		this.localAlignment=localAlignment;

		int readOneLength = localAlignment.getOriginalSequences().get(0).getLength();
		int readTwoLength = localAlignment.getOriginalSequences().get(1).getLength();
		int alignmentBlockLength = localAlignment.getLength();
		// rare case if there's no alignment block at all.
		if (alignmentBlockLength==0) {
			int maxLen = Math.max(readOneLength, readTwoLength);
			int [] [] result = new int [2][maxLen];
			for (int i=0; i<readOneLength; i++)
				result[0][i]=i+1;
			for (int i=0; i<readTwoLength; i++)
				result[1][i]=i+1;
			this.index=result;
			return;
		}

		int queryStart = localAlignment.getIndexInQueryAt(1);
		int targetStart = localAlignment.getIndexInTargetAt(1);
		int queryEnd = localAlignment.getIndexInQueryAt(localAlignment.getLength());
		int targetEnd = localAlignment.getIndexInTargetAt(localAlignment.getLength());
		int trailingBasesOne =   readOneLength - queryEnd;
		int trailingBasesTwo =   readTwoLength - targetEnd;

		// the end position difference between read one and read two - positive means read one is "ahead" of read two.
		// int postAlignmentBlockOffSet = queryEnd - targetEnd;

		int length1 = Math.max(queryStart-1, targetStart-1);

		int length3 = Math.max(trailingBasesOne, trailingBasesTwo);
		int consensusLength = length1+alignmentBlockLength+length3;
		int [] [] index = new int [2][consensusLength];

		index=getPreAlignmentIndexes(queryStart, targetStart, index);
		index= getAlignmentIndexes(queryStart, targetStart, index, localAlignment);
		index=getPostAlignmentIndexes(queryEnd, targetEnd, trailingBasesOne, trailingBasesTwo, localAlignment, index, readOneLength, readTwoLength);
		this.index=index;
	}

	@Override
	public String toString() {
		return Arrays.toString(index[0]) + "\n" + Arrays.toString(index[1]);
	}

	private int [] []  getPreAlignmentIndexes (final int queryStart, final int targetStart, final int [] [] index) {
		// the start position difference between read one and read two - positive means read one has more bases, negative means read two has more bases.
		// int preAlignmentBlockOffSet = queryStart - targetStart;

		// 0 based.
		int maxLength=Math.max(queryStart, targetStart);
		for (int i=0; i<maxLength-1; i++) {
			int ts = targetStart- maxLength+i;
			int qs = queryStart - maxLength+i;
			if (ts<0) ts=-1;
			if (qs<0) qs=-1;
			index[0][i]=qs+1;// make one based.
			index[1][i]=ts+1;
		}
		return index;
	}

	private int [] [] getAlignmentIndexes (final int queryStart, final int targetStart, final int [] [] index, final SequencePair<DNASequence, NucleotideCompound> localAlignment) {
		int maxLength=Math.max(queryStart, targetStart);  // where to start filling in the local alignment
		int alignmentLength = localAlignment.getLength();

		AlignedSequence<DNASequence,NucleotideCompound>  alignedSeq1 = localAlignment.getAlignedSequence(1);
		AlignedSequence<DNASequence,NucleotideCompound>  alignedSeq2 = localAlignment.getAlignedSequence(2);

		for (int i=0; i<alignmentLength; i++) {
			int idx = i+maxLength-1;
			int b1 = localAlignment.getIndexInQueryAt(i+1); // one based index
			int b2 = localAlignment.getIndexInTargetAt(i+1);

			if (alignedSeq1.isGap(i+1)) b1=0;
			if (alignedSeq2.isGap(i+1)) b2=0;
			index[0][idx]=b1;
			index[1][idx]=b2;
		}
		return index;
	}

	/*
	private int [] [] getPostAlignmentIndexes (final int queryEnd, final int targetEnd, final int trailingBasesOne, final int trailingBasesTwo,  final SequencePair<DNASequence, NucleotideCompound> localAlignment, final int [] [] index) {
		int queryStart = localAlignment.getIndexInQueryAt(localAlignment.getLength());
		int targetStart = localAlignment.getIndexInTargetAt(localAlignment.getLength());

		int maxLength=Math.max(queryStart, targetStart);

		// fill in each individually
		for (int i=0; i<trailingBasesOne; i++) {
			int queryIndex = maxLength+i;
			index[0][queryIndex-1]=queryStart+i;
		}
		// target.
		for (int i=0; i<trailingBasesTwo; i++) {
			int targetIndex = maxLength+i;
			index[1][targetIndex-1]=targetStart+i;
		}

		return index;
	}
	*/

	private int [] [] getPostAlignmentIndexes (final int queryEnd, final int targetEnd, final int trailingBasesOne, final int trailingBasesTwo,  final SequencePair<DNASequence, NucleotideCompound> localAlignment, final int [] [] index,
			final int readOneLength, final int readTwoLength) {
		// iterate backwards?
		int consensusSize=index[0].length;
		int r1ConsensusPosition = consensusSize-trailingBasesOne;
		int r2ConsensusPosition = consensusSize -trailingBasesTwo;
		int startConsensusPosition = Math.min(r1ConsensusPosition, r2ConsensusPosition);
		int r1FirstBase = readOneLength - trailingBasesOne+1;
		int r2FirstBase = readTwoLength - trailingBasesTwo+1;

		for (int i=0; i<trailingBasesOne; i++)
			index[0][startConsensusPosition+i]=i+r1FirstBase;

		for (int i=0; i<trailingBasesTwo; i++)
			index[1][startConsensusPosition+i]=i+r2FirstBase;


		return index;
	}


	/**
	 * Get the position in read one, based on the read consensus position.  This is one-based.
	 * @param consensusPosition
	 * @return
	 */
	public int getReadOnePosition(final int consensusPosition) {
		return this.index[0][consensusPosition-1];
	}

	public int getReadTwoPosition (final int consensusPosition) {
		return this.index[1][consensusPosition-1];
	}

	/**
	 * Get the position in read one (1) or two (2) based on the consensus position.
	 * @param readNum
	 * @param consensusPosition
	 * @return
	 */
	public int getReadPosition (final int readNum, final int consensusPosition) {
		if (readNum<1 || readNum >2) throw new IllegalArgumentException("readNum should be 1 or 2");
		return this.index[readNum-1][consensusPosition-1];
	}

	public int getConsensusLength () {
		return this.index[0].length;
	}

}
