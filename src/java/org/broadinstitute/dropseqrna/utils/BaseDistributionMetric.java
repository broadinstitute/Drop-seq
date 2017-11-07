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
