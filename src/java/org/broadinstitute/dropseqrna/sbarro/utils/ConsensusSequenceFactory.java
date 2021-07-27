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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqRecord;
import htsjdk.samtools.util.SequenceUtil;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.SubstitutionMatrixHelper;
import org.biojava.nbio.alignment.template.PairwiseSequenceAligner;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;

public class ConsensusSequenceFactory {

	private static ConsensusSequenceFactory instance = null;
	private final SubstitutionMatrix<NucleotideCompound> matrix;
	private final SimpleGapPenalty gapP;

	protected ConsensusSequenceFactory () {
		this.matrix = SubstitutionMatrixHelper.getNuc4_4();
		SimpleGapPenalty gapPenality = new SimpleGapPenalty();
		gapPenality.setOpenPenalty(5);
		gapPenality.setExtensionPenalty(2);
		this.gapP=gapPenality;
	}

	public static ConsensusSequenceFactory getInstance() {
		if (instance==null)
			synchronized (ConsensusSequenceFactory.class) {
				if (instance==null)
					instance = new ConsensusSequenceFactory();
			}
		return instance;
	}

	/**
	 * For 2 DNA sequences, align them and generate a consensus sequence result.
	 * This attempts both consensus of both the two sequences and the reverse compliment of the second sequence to find the best match.
	 * The two alignments are compared by score, with the better score becoming the returned alignment.
	 * @param seq1
	 * @param seq2
	 * @return
	 */
	public ConsensusSequence getConsensusSequence (final String seq1, final String seq2, final boolean assumeRC) {
		ConsensusSequence alignedPair = getBestAlignedPair(seq1, seq2, assumeRC);
		return alignedPair;
	}

	public ConsensusSequence getConsensusSequence (final FastqRecord seq1, final FastqRecord seq2, final boolean assumeRC) {
		String n1 = SequenceUtil.getSamReadNameFromFastqHeader(seq1.getReadName());
		String n2 = SequenceUtil.getSamReadNameFromFastqHeader(seq2.getReadName());
		if (!n1.equals(n2))
			throw new IllegalStateException("Paired reads don't have the same  name: " + n1 + " " + n2);
		ConsensusSequence alignedPair = getBestAlignedPair(seq1.getReadString(), seq2.getReadString(), assumeRC);
		alignedPair.addReadBaseQualities(seq1.getBaseQualityString(), seq2.getBaseQualityString());
		return alignedPair;
	}

	public ConsensusSequence getConsensusSequence (final SAMRecord r1, final SAMRecord r2, final boolean assumeRC) {
		if (!r1.getReadName().equals(r2.getReadName()))
				throw new IllegalStateException("Paired reads don't have the same  name: " + r1.getReadName() + " " + r2.getReadName());
		ConsensusSequence alignedPair = getBestAlignedPair(r1.getReadString(), r1.getReadString(),assumeRC);
		alignedPair.addReadBaseQualities(r1.getBaseQualityString(), r2.getBaseQualityString());
		return alignedPair;
	}

	/**
	 * Performs global alignment to both the default sequence and the reverse compliment of the second sequence, and selects the best result by score.
	 *
	 * @param seq1 The first read sequence
	 * @param seq2 The second read sequence
	 * @return The parsimonious pair of aligned sequences.
	 */
	/*
	private ConsensusSequence getBestAlignedPair (final String seq1, final String seq2) {
		PairwiseSequenceAligner <DNASequence,NucleotideCompound>  defaultAlignment = Alignments.getPairwiseAligner(AlignmentUtils.getDNASequenceFromString(seq1),
				AlignmentUtils.getDNASequenceFromString(seq2), PairwiseSequenceAlignerType.GLOBAL, gapP, matrix);

		PairwiseSequenceAligner <DNASequence,NucleotideCompound> reverseComplimentAlignment = Alignments.getPairwiseAligner(AlignmentUtils.getDNASequenceFromString(seq1),
				AlignmentUtils.getDNASequenceFromString(SequenceUtil.reverseComplement(seq2)), PairwiseSequenceAlignerType.GLOBAL, gapP, matrix);

		double score1 = defaultAlignment.getScore();
		double score2 = reverseComplimentAlignment.getScore();
		// pick the best score and set up the consensus object.
		ConsensusSequence result = null;

		if (score1 > score2) {
			PairwiseSequenceAligner <DNASequence,NucleotideCompound>  localAlignment = Alignments.getPairwiseAligner(AlignmentUtils.getDNASequenceFromString(seq1),
					AlignmentUtils.getDNASequenceFromString(seq2), PairwiseSequenceAlignerType.LOCAL, gapP, matrix);
			result = new ConsensusSequence(defaultAlignment.getPair(),localAlignment.getPair(), score1, false);
		}
		else {
			PairwiseSequenceAligner <DNASequence,NucleotideCompound>  localAlignment = Alignments.getPairwiseAligner(AlignmentUtils.getDNASequenceFromString(seq1),
					AlignmentUtils.getDNASequenceFromString(SequenceUtil.reverseComplement(seq2)), PairwiseSequenceAlignerType.LOCAL, gapP, matrix);
			result = new ConsensusSequence(reverseComplimentAlignment.getPair(), localAlignment.getPair(), score2, true);
		}
		return result;

	}
	*/

	/**
	 * If assumeRC is true, then the 2nd read will be reverse complimented before consensus is generated.
	 * Otherwise, the 2nd read sequence and RC will be tested, resulting in 2x as much computational work.
	 * @param seq1
	 * @param seq2
	 * @param assumeRC
	 * @return
	 */
	private ConsensusSequence getBestAlignedPair (final String seq1, final String seq2, final boolean assumeRC) {

		if (seq1.length()==0 || seq2.length()==0) throw new IllegalArgumentException("At least one sequence was length 0.  Can not generated consensus.  R1 [" + seq1 +"] R2 ["+ seq2 +"]");

		PairwiseSequenceAligner <DNASequence,NucleotideCompound> reverseComplimentAlignment = Alignments.getPairwiseAligner(AlignmentUtils.getDNASequenceFromString(seq1),
				AlignmentUtils.getDNASequenceFromString(SequenceUtil.reverseComplement(seq2)), PairwiseSequenceAlignerType.LOCAL, gapP, matrix);

		if (assumeRC) {
			ConsensusSequence result = new ConsensusSequence(reverseComplimentAlignment.getPair(), 0, true);
			return (result);
		}

		// if you're not assuming, then do both alignments and get the best score.
		PairwiseSequenceAligner <DNASequence,NucleotideCompound>  defaultAlignment = Alignments.getPairwiseAligner(AlignmentUtils.getDNASequenceFromString(seq1),
				AlignmentUtils.getDNASequenceFromString(seq2), PairwiseSequenceAlignerType.LOCAL, gapP, matrix);


		double score1 = defaultAlignment.getScore();
		double score2 = reverseComplimentAlignment.getScore();
		// pick the best score and set up the consensus object.
		ConsensusSequence result = null;

		if (score1 > score2)
			result = new ConsensusSequence(defaultAlignment.getPair(), score1, false);
		else
			result = new ConsensusSequence(reverseComplimentAlignment.getPair(), score2, true);
		return result;

	}

}
