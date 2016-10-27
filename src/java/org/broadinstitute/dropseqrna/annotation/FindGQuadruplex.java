package org.broadinstitute.dropseqrna.annotation;

import java.util.regex.Pattern;

/**
 * Given a sequence, find a list of GQuadruplex intervals.
 * d(G3+ N1-7 G3+ N1-7 G3+ N1-7 G3+), where N is any base (including guanine)
 * @author nemesh
 *
 */

public class FindGQuadruplex {

	private static Pattern pattern;

	public FindGQuadruplex () {

	}

	/**
	 * This finds all non-overlapping G-Quadraplex substrings.
	 * If there are ambiguous positions for the G-Quadraplex, the first one is recovered.
	 * @param seqName
	 * @param seq
	 */

}
