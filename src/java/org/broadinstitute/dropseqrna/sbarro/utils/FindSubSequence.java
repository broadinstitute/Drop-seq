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

import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.alignment.SubstitutionMatrixHelper;
import org.biojava.nbio.alignment.template.SequencePair;
import org.biojava.nbio.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;

/**
 * Given a read sequence and a target sequence, find the offset of the target sequence in the read sequence.
 * For example, given a sequence "AAGGTTTTGAAAAA", find "TTT"'s offset in this sequence.
 * @author nemesh
 *
 */
public class FindSubSequence {


	private final DNASequence target;
	private final SubstitutionMatrix<NucleotideCompound> matrix;
	private final SimpleGapPenalty gapP;

	/**
	 * Construct a global alignment tool that looks for a target sequence in submitted query sequences.
	 * @param targetSequence The sequence to search for in queries.
	 * @param gapOpenPenalty The cost for opening a gap
	 * @param gapExtensionPenality The cost for extending a gap
	 */
	public FindSubSequence (final String targetSequence, final int gapOpenPenalty, final int gapExtensionPenality) {
		this.target=AlignmentUtils.getDNASequenceFromString(targetSequence);
		this.matrix = SubstitutionMatrixHelper.getNuc4_4();
		SimpleGapPenalty gapPenality = new SimpleGapPenalty();
		gapPenality.setOpenPenalty(gapOpenPenalty);
		gapPenality.setExtensionPenalty(gapExtensionPenality);
		this.gapP=gapPenality;
	}

	/**
	 * Construct a global alignment tool that looks for a target sequence in submitted query sequences.
	 * Initialized with sensible gap penalty (5) and extension penalty (2).
	 * @param targetSequence The sequence to search for in queries.
     */
	public FindSubSequence (final String targetSequence) {
		this(targetSequence, 5, 2);
	}

	/**
	 * Given a query sequence, find the location in the target using a global alignment strategy
	 * @param query
	 * @return
	 */
	public SubSequenceResultI findSequenceGlobalAlignment (final String query) {
		DNASequence q = AlignmentUtils.getDNASequenceFromString(query);
		SequencePair<DNASequence,NucleotideCompound> p = Alignments.getPairwiseAlignment(q, target,PairwiseSequenceAlignerType.GLOBAL, gapP, matrix);
		return new SubSequenceResultGlobalAlignment(p);
	}

	/**
	 * Given a query sequence, find the location in the target using a global alignment strategy
	 * @param query
	 * @return
	 */
	public SubSequenceResultI findSequenceLocalAlignment (final String query) {
		DNASequence q = AlignmentUtils.getDNASequenceFromString(query);
		SequencePair<DNASequence,NucleotideCompound> p2 = Alignments.getPairwiseAlignment(q, target,PairwiseSequenceAlignerType.LOCAL, gapP, matrix);
		return new SubSequenceResultLocalAlignment(p2);
	}






}
