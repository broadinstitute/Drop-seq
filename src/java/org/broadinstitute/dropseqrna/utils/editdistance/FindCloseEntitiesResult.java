package org.broadinstitute.dropseqrna.utils.editdistance;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class FindCloseEntitiesResult <T,M> {
	private final Map<T, List<T>> entityMap;
	private final Set<M> collapseMetric;

	public FindCloseEntitiesResult () {
		this.entityMap=new HashMap<T, List<T>>();
		this.collapseMetric=new HashSet<M>();
	}

	public void addMapping (Map<T, List<T>> entityMap) {
		for (T key: entityMap.keySet())
			addMapping(key, entityMap.get(key));
	}
			
	public void addMapping (T key, Collection<T> value) {
		entityMap.put(key, new ArrayList<T>(value));
	}
	
	public void addMetrics (Collection<M> metrics) {
		this.collapseMetric.addAll(metrics);
	}
	
	public void addMetrics (M metrics) {
		this.collapseMetric.add(metrics);
	}
	
	public Map<T, List<T>> getEntityMap() {
		return entityMap;
	}

	public Set<M> getCollapseMetric() {
		return collapseMetric;
	}
	
	
	
}
