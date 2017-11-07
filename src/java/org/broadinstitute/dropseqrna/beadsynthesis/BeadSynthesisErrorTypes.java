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
package org.broadinstitute.dropseqrna.beadsynthesis;

public enum BeadSynthesisErrorTypes {
	NO_ERROR("NO_ERROR"),
	SYNTH_MISSING_BASE("SYNTH_MISSING_BASE"),
	SINGLE_UMI("SINGLE_UMI"),
	PRIMER("PRIMER"),
	FIXED_FIRST_BASE("FIXED_FIRST_BASE"),
	OTHER_ERROR("OTHER_ERROR");
	
	
	private final String text;
	
	private BeadSynthesisErrorTypes(String text) {
		this.text=text;
	}
	
	public String toString () {
		return this.text;
	}
}
