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

	int getCount(Character base) {
		return (map.get(base));
	}

	public String toString() {
		return this.map.toString();
	}

}
