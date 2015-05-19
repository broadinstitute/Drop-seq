package org.broadinstitute.dropseqrna.readtrimming;

/**
 * Finds runs of PolyAs (lots of A's in a row) in a sequence.
 * Identifies the first base of the run that is polyA.
 * @author nemesh
 *
 */
public class PolyAFinder {

	private int minNumBases;
	private int numMismatchBases;
	
	public PolyAFinder (int minNumBases, int numMismatchBases) {
		this.minNumBases = minNumBases;
		this.numMismatchBases = numMismatchBases;
	}
	
	/**
	 * Gets the first base of a polyA run in the sequence, or -1 if there is no polyA run of at least minNumBases.
	 * @param sequence the sequence to test.
	 * @return The index of the first base of the polyA run.  This is a 0 based index!
	 */
	public int getPolyAStart(String sequence) {
		sequence=sequence.toUpperCase();
		char [] seq = sequence.toCharArray();
		
		boolean inRun=false;
		int numInRun=0;
		int numError=0;
		int runStartPos=-1;
		int bestRunStart=-1;
		int bestNumInRun=-1;
		for (int i=0; i<seq.length; i++)  {
			char s = seq[i];
			// keep the run going if you're in one.
			if (s=='A' && inRun) {
				numInRun++;
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
					}
					runStartPos=-1; // you have too many mismatches, end the run and reset.
					numInRun=0;
					numError=0;
					inRun=false;
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
		}
		
		// if you have enough bases to consider this a run, great.
		if (bestNumInRun>=this.minNumBases) return bestRunStart;
		// otherwise, return -1 because you don't have a start postion for a run.
		return -1;
	}
	
	
	
}
