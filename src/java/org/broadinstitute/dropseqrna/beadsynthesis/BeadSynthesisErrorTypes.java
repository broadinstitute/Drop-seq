package org.broadinstitute.dropseqrna.beadsynthesis;

public enum BeadSynthesisErrorTypes {
	NO_ERROR("NO_ERROR"),
	POLY_T_ERROR("POLY_T_ERROR"),
	SINGLE_UMI("SINGLE_UMI"),
	OTHER_ERROR("OTHER_ERROR");
	
	private final String text;
	
	private BeadSynthesisErrorTypes(String text) {
		this.text=text;
	}
	
	public String toString () {
		return this.text;
	}
}
