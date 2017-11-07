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
