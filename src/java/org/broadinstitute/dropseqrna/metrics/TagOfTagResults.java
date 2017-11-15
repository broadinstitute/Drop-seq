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
