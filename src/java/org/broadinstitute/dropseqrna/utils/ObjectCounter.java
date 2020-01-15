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

import org.broadinstitute.dropseqrna.TranscriptomeException;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Keeps track of the number of times an object has been seen.
 * @author nemesh
 *
 * @param <T> The type of object to count.  Object needs to have an equals method implemented!
 */
public class ObjectCounter<T extends Comparable<T>> {

	private Map<T, Integer> countMap;

	public ObjectCounter () {
		countMap = new HashMap<>();
	}

	/**
	 * Make a shallow copy of the input object counter.
	 */
	public ObjectCounter (final ObjectCounter<T> counter) {
		this();
		for (T key: counter.getKeys())
			this.countMap.put(key, counter.getCountForKey(key));
	}

	public boolean hasKey (final T object) {
		return this.countMap.containsKey(object);
	}

	public void increment (final T object) {
		incrementByCount(object, 1);
	}

	public void decrement (final T object) {
		incrementByCount(object, -1);
	}

	public void decrementByCount (final T object, final int count) {
		int count2=count *-1;
		incrementByCount(object, count2);
	}

	/**
	 * Add the contents of an ObjectCounter of the same type to this object.
	 * @param object
	 */
	public void increment (final ObjectCounter<T> object) {
		for (T key : object.getKeys()) {
			Integer count = object.getCountForKey(key);
			this.incrementByCount(key, count);
		}
	}

	public void clear() {
		this.countMap.clear();
	}

	public void incrementByCount (final T object, final int size) {
		Integer count = countMap.get(object);
		if (count==null)
			countMap.put(object, size);
		else {
			count+=size;
			countMap.put(object, count);
		}
	}

	public void setCount(final T object, final int count) {
		countMap.put(object, count);
	}

	public void remove (final T object) {
		this.countMap.remove(object);
	}

	public Collection<T> getKeys () {
		return countMap.keySet();
	}

	public int getSize () {
		return countMap.size();
	}

	public int getCountForKey (final T key) {
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

	public int getNumberOfSize(final int size) {
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

	/**
	 * Get the object with the fewest counts.
	 * @return
	 */
	public T getMin() {
		T min=null;
		int minCount=Integer.MAX_VALUE;
		for (T key: this.getKeys()) {
			int count=this.getCountForKey(key);
			if (count<minCount) {
				min = key;
				minCount=count;
			}
		}
		return (min);

	}

	public List<T> getKeysOrderedByCount (final boolean decreasing) {
		Map<Integer, List<T>> reversed= this.getReverseMapping();
		List<Integer> counts = new ArrayList<>(reversed.keySet());
		Collections.sort(counts);
		if (decreasing) Collections.reverse(counts);
		List<T> keys = new ArrayList<>();
		for (int i: counts) {
			List<T> t = reversed.get(i);
			Collections.sort(t);
			keys.addAll(t);
		}
		return (keys);
	}

	public Map<Integer, List<T>> getReverseMapping () {
		Map<Integer, List<T>> result = new HashMap<>(this.countMap.size());
		for (T key : this.countMap.keySet()) {
			Integer v =getCountForKey(key);
			List<T> l = result.get(v);
			if (l==null)
				l=new ArrayList<>();
			l.add(key);
			result.put(v, l);
		}
		return result;
	}

	/**
	 * Filters this counter to that only entries with at least <count> number of reads remain.
	 */
	public void filterByMinCount (final int count) {
		Map<T, Integer> result = new HashMap<>();
		for (T key: this.countMap.keySet()) {
			Integer value = countMap.get(key);
			if (value>=count)
				result.put(key, value);
		}
		this.countMap=result;
	}

	/**
	 * Subset this set of counts to a subset of the keys.
	 * @param keys A collection of keys to restrict the data to.
	 */
	public void subset (final Set<T> keys) {
		Map<T, Integer> newMap = new HashMap<>();
		for (T k: countMap.keySet())
			if (keys.contains(k))
				newMap.put(k, this.countMap.get(k));

		this.countMap=newMap;
	}

	@Override
	public String toString () {
		return this.countMap.toString();
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((countMap == null) ? 0 : countMap.hashCode());
		return result;
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		ObjectCounter other = (ObjectCounter) obj;
		if (countMap == null) {
			if (other.countMap != null)
				return false;
		} else if (!countMap.equals(other.countMap))
			return false;
		return true;
	}

    @SuppressWarnings("unchecked")
    public boolean countMapsEqual(final ObjectCounter<T> otherCounter) {
        if (this.equals(otherCounter))
            return true;
        if (otherCounter == null)
            return false;

        List<T> keys = getKeysOrderedByCount(true);
        List<T> otherKeys = otherCounter.getKeysOrderedByCount(true);
        if (otherKeys.size() != keys.size())
            return false;

        for (int idx=0; idx<keys.size(); idx++) {
            T key = keys.get(idx);
            T otherKey = otherKeys.get(idx);
            if (!otherKey.equals(key) || otherCounter.getCountForKey(otherKey) != this.getCountForKey(key))
                return false;
        }

        return true;
    }

    /**
     * @param reportFile The report file to read.
     */
    public static ObjectCounter<String> readReportFile(final File reportFile) {
        ObjectCounter <String> counter = new ObjectCounter<>();

        try {
            BufferedReader input = new BufferedReader(new FileReader(reportFile));
            try {
                String line;
                while ((line = input.readLine()) != null) {
                    line = line.trim();
                    if (line.startsWith("#")) {
                        continue;
                    }
                    String[] strLine = line.split("\t");
                    int count = Integer.parseInt(strLine[0]);
                    String key = strLine[1];
                    counter.incrementByCount(key, count);
                }
            } finally {
                input.close();
            }
        } catch (IOException ex) {
            throw new TranscriptomeException("Error reading the file: " + reportFile.getAbsolutePath());
        }

        return counter;
    }
}
