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
package org.broadinstitute.dropseqrna.metrics;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class TagOfTagResults<KEY,VALUE> {

	private Map<KEY, Set<VALUE>> result = null;
	
	// values will be used many times in the result map, so use references to those values.
	private Map<VALUE, VALUE> valueCache = null;
	
	public TagOfTagResults () {
		result = new HashMap<KEY, Set<VALUE>>();
		valueCache = new HashMap<VALUE, VALUE>();
	}
	
	public Set<KEY> getKeys () {
		return result.keySet();
	}
	
	public Set<VALUE> getValues (String key) {
		return result.get(key);
	}
	
	public Integer getCount(String key) {
		return result.get(key).size();
	}
		
	public void addEntry (KEY key, VALUE value) {
		Set<VALUE> values = result.get(key);
		if (values==null) {
			values = new HashSet<VALUE>();
			result.put(key, values);
		}
		VALUE v = checkCache(value);
		values.add(v);
	}
	
	public void addEntries(KEY key, Collection<VALUE> value) {
		for (VALUE v: value) {
			addEntry(key, v);
		}
	}
	
	private VALUE checkCache(VALUE key) {
		VALUE v = this.valueCache.get(key);
		if (v!=null) return (v);
		// otherwise, add and return.
		this.valueCache.put(key, key);
		return (checkCache(key));
	}
	
}
