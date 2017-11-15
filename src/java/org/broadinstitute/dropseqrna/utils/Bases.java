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
