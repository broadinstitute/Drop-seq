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

import htsjdk.samtools.metrics.MetricBase;

import java.util.HashMap;
import java.util.Map;

public class BaseDistributionMetric extends MetricBase {

	private Map<Character, Integer> map = null;

	public BaseDistributionMetric() {
		map = new HashMap<Character, Integer>();
		map.put(Bases.A.getBase(), 0);
		map.put(Bases.C.getBase(), 0);
		map.put(Bases.G.getBase(), 0);
		map.put(Bases.T.getBase(), 0);
		map.put(Bases.N.getBase(), 0);
	}

	void addBase(Character base) {
		int count = map.get(base);
		count++;
		map.put(base, count);
	}

	public int getCount(Character base) {
		return (map.get(base));
	}
	
	public int getTotalCount () {
		int count=0;
		for (Bases b: Bases.values()) {
			char bb = b.getBase();
			count+=map.get(bb);
		}
		return count;
	}

	public String toString() {
		return this.map.toString();
	}

}
