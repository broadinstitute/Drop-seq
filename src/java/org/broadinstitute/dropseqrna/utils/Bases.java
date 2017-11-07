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
package org.broadinstitute.dropseqrna.utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public enum Bases {
	A('A'),
	C('C'),
	G('G'),
	T('T'),
	N('N');

	private Character base;

	Bases(final Character base) {
		this.base=base;
	}

	public Character getBase() {
		return this.base;
	}

	/**
	 * Get one of the bases that isn't this one.
	 * @param b The base to not include as a possible outcome
	 * @param excludeN Should the base N be excluded as a possible draw?
	 * @return
	 */
	public Bases getOtherBase(final Bases b, final boolean excludeN) {
		List<Bases> all = new ArrayList<Bases>(Arrays.asList(Bases.values()));
		if (excludeN)
			all.remove(Bases.N);
		all.remove(b);
		Collections.shuffle(all);
		return all.get(0);
	}

}
