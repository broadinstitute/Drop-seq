package org.broadinstitute.dropseqrna.utils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class BaseDistributionMetricCollection {
	
	private Map<Integer, BaseDistributionMetric> collection = null;
	
	public BaseDistributionMetricCollection() {
		collection = new HashMap<Integer, BaseDistributionMetric>();
	}
	
	public void addBase (char base, int position) {
		BaseDistributionMetric m = this.collection.get(position);
		if (m==null) {
			m = new BaseDistributionMetric();
			this.collection.put(position, m);
		}
		m.addBase(base);
	}
	
	public void addBases (String bases) {
		char [] b = bases.toCharArray();
		for (int i=0; i<b.length; i++) {
			addBase(b[i], i);
		}
	}
	
	public void addBases (char [] bases) {
		for (int i=0; i<bases.length; i++) {
			addBase(bases[i], i);
		}
	}
	
	public void addBases (byte [] bases) {
		for (int i=0; i<bases.length; i++) {
			char b = (char) bases[i];
			addBase(b, i);
		}
	}
	
	public BaseDistributionMetric getDistributionAtPosition (int position) {
		return this.collection.get(position);
	}
	
	public List<Integer> getPositions () {
		List<Integer> r = new ArrayList<Integer>(this.collection.keySet());
		Collections.sort(r);
		return (r);
	}
}
