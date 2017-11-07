/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.dropseqrna.readtrimming;

import htsjdk.samtools.SAMRecord;

/**
 * Finds runs of PolyAs (lots of A's in a row) in a sequence.
 * Identifies the first base of the run that is polyA.
 * @author nemesh
 *
 */
public class SimplePolyAFinder implements PolyAFinder {

	private int minNumBases;
	private int numMismatchBases;

	public SimplePolyAFinder(int minNumBases, int numMismatchBases) {
		this.minNumBases = minNumBases;
		this.numMismatchBases = numMismatchBases;
	}

    /**
	 * Gets the first base of a polyA run in the sequence, or -1 if there is no polyA run of at least minNumBases.
	 * @param rec the read to test
	 * @return The index of the first base of the polyA run.  This is a 0 based index!
	 */
	@Override
    public PolyARun getPolyAStart(final SAMRecord rec) {
        return getPolyAStart(rec.getReadString());
    }

    public PolyARun getPolyAStart(final String readString) {
        char [] seq = readString.toUpperCase().toCharArray();

		boolean inRun=false;
		int numInRun=0;
		int numError=0;
		int runStartPos=-1;
		int bestRunStart=-1;
		int bestNumInRun=-1;
        int bestRunLength= -1;
        int nextRunStartPos = -1;
		for (int i=0; i<seq.length; i++)  {
			char s = seq[i];
			// keep the run going if you're in one.
			if (s=='A' && inRun) {
				numInRun++;
                if (seq[i-1] != 'A' && nextRunStartPos == -1) {
                    // If hit a mismatch but continuing, remember where the first A was after the first mismatch
                    nextRunStartPos = i;
                }
			}
		
			// Get the run going if you're starting one.  
			if (s=='A' & !inRun) {
				numInRun++;
				inRun=true;
				runStartPos=i;
				
			} 
			// if you're in a run and this isn't an A, then you accumulate errors, and maybe leave the run if you have too many errors.
			// if you haven't accumulated too many errors, keep the run going as if it was an OK base.
			if (s!='A' & inRun) {
				numError++;
				numInRun++;
				if (numError>this.numMismatchBases) {
					if ((numInRun-numError)>bestNumInRun) {
						bestRunStart=runStartPos;
						bestNumInRun=numInRun-numError;
                        bestRunLength= i - runStartPos;
					}
					runStartPos=-1; // you have too many mismatches, end the run and reset.
					numInRun=0;
					numError=0;
					inRun=false;
                    // Jump back to just before the the first A after the first mismatch run in the run being cleared.
                    if (nextRunStartPos != -1) {
                        i = nextRunStartPos - 1;
                        nextRunStartPos = -1;
                    }
				}
				
				
			}
			if (s!='A' & !inRun) {
				// you really don't do much here.
			}
				

		}
		
		// finally, if you're in a run, see if that was better than the previous run when you ran out of sequence.
		if (inRun && (numInRun-numError)>bestNumInRun) {
			bestRunStart=runStartPos;
			bestNumInRun=numInRun-numError;
            bestRunLength = seq.length - runStartPos;
		}
		
		// if you have enough bases to consider this a run, great.
		if (bestNumInRun>=this.minNumBases) return new PolyARun(bestRunStart, bestRunLength);
		// otherwise, return -1 because you don't have a start postion for a run.
		return PolyARun.NO_MATCH_RUN;
	}
	
	
	
}
