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
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Keeps track of the number of times an object has been seen.
 * @author nemesh
 *
 * @param <T> The type of object to count.  Object needs to have an equals method implemented!
 */
public class ObjectCounter<T extends Comparable<T>> {

	private Map<T, Integer> countMap;
	
	public ObjectCounter () {
		countMap = new HashMap<T, Integer>();
	}
	
	/**
	 * Make a shallow copy of the input object counter.
	 */
	public ObjectCounter (ObjectCounter<T> counter) {
		this();
		for (T key: counter.getKeys()) {
			this.countMap.put(key, counter.getCountForKey(key));
		}
	}
	
	public boolean hasKey (T object) {
		return this.countMap.containsKey(object);
	}
	
	public void increment (T object) {
		incrementByCount(object, 1);
	}
	
	/**
	 * Add the contents of an ObjectCounter of the same type to this object. 
	 * @param object
	 */
	public void increment (ObjectCounter<T> object) {
		for (T key : object.getKeys()) {
			Integer count = object.getCountForKey(key);
			this.incrementByCount(key, count);
		}
	}
	
	public void clear() {
		this.countMap.clear();
	}
	
	public void incrementByCount (T object, int size) {
		Integer count = countMap.get(object);
		if (count==null) {
			countMap.put(object, size);
		} else {
			count+=size;
			countMap.put(object, count);
		}
	}
	
	public void setCount(T object, int count) {
		countMap.put(object, count);
	}
	
	public void remove (T object) {
		this.countMap.remove(object);
	}
	
	public Collection<T> getKeys () {
		return countMap.keySet();
	}
	
	public int getSize () {
		return countMap.size();
	}
	
	public int getCountForKey (T key) {
		Integer count = countMap.get(key);
		if (count==null) return 0;
		return count;
	}
	
	public Collection<Integer> getCounts() {
		return (this.countMap.values());
	}
	
	public int getTotalCount() {
		int result = 0;
		for (T key: this.countMap.keySet()) {
			int t = getCountForKey(key);
			result+=t;
		}
		return result;
	}
	
	public int getNumberOfSize(int size) {
		int result = 0;
		for (T key: this.countMap.keySet()) {
			int t = getCountForKey(key);
			if (t==size) result++;
		}
		return result;
	}
	
	public T getMode () {
		T max = null;
		int maxCount=0;
		for (T key: this.getKeys()) {
			int count=this.getCountForKey(key);
			if (count>maxCount) {
				max = key;
				maxCount=count;
			}
		}
		return (max);
	}
	// NOTE: 
	// for keys with the same number of items, object ordering is undefined.  Need to fix this to break ties by the T's natural ordering.	
	public List<T> getKeysOrderedByCount (boolean decreasing) {
		Map<Integer, List<T>> reversed= this.getReverseMapping();
		List<Integer> counts = new ArrayList<Integer>(reversed.keySet());
		Collections.sort(counts);
		if (decreasing) Collections.reverse(counts);
		List<T> keys = new ArrayList<T>();
		for (int i: counts) {
			List<T> t = reversed.get(i);
			Collections.sort(t);
			keys.addAll(t);
		}
		return (keys);
	}
	
	public Map<Integer, List<T>> getReverseMapping () {
		Map<Integer, List<T>> result = new HashMap<Integer, List<T>>(this.countMap.size());
		for (T key : this.countMap.keySet()) {
			Integer v =getCountForKey(key);
			List<T> l = result.get(v);
			if (l==null) {
				l=new ArrayList<T>();
			} 
			l.add(key);
			result.put(v, l);
		}
		return result;
	}
	
	/**
	 * Filters this counter to that only entries with at least <count> number of reads remain. 
	 */
	public void filterByMinCount (int count) {
		Map<T, Integer> result = new HashMap<T, Integer>();
		for (T key: this.countMap.keySet()) {
			Integer value = countMap.get(key);
			if (value>=count) {
				result.put(key, value);
			}
		}
		this.countMap=result;
	}
	
	// public boolean equals()
	
	public String toString () {
		return this.countMap.toString();
	}
	
	
	
	
}
