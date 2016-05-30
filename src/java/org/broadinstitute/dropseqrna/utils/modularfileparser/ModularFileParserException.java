package org.broadinstitute.dropseqrna.utils.modularfileparser;

public class ModularFileParserException extends RuntimeException {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1605829591524594060L;

	public ModularFileParserException(final String s) {
		super(s);
	}

	public ModularFileParserException(final String s, final Throwable throwable) {
		super(s, throwable);
	}
}
