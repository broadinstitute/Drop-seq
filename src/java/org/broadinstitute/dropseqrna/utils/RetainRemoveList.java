/*
 * MIT License
 *
 * Copyright 2018 Broad Institute
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
import java.util.List;
import java.util.Set;

/**
 * Utility class to retain/remove elements from a list.
 * @author nemesh
 *
 * @param <T>
 */
public class RetainRemoveList <T> {

	/**
	 * Given a list of original elements and an optional set to remove or retain, return elements from the original list retaining their order, removing elements
	 * that are specified to be removed, and retaining elements that are specified to be retained.
	 * Retain takes priority over remove if the element is in both sets.
	 * If the element is in neither set, then the result depends on if the retain and remove sets are populated.
	 * If retain is populated and the element isn't in it, do not retain this element
	 * If retain is size 0 and the element isn't in the remove set, retain this element
	 * @param elements The ordered list of element to test
	 * @param retain A set of elements to retain
	 * @param remove A set of elements to remove.
	 * @return The ordered list of elements to return
	 */
	public List<T> getElementsToRetain (final List<T> elements, final Set<T> toRemove, final Set<T> toRetain) {
		List<T> result = new ArrayList<>();

		for (T h: elements)
			if (retainElement(h, toRetain, toRemove))
				result.add(h);
		return (result);
	}

	/**
	 * Should this element be retained, based on the retain and remove sets?
	 * Retain takes priority over remove if the element is in both sets.
	 * If the element is in neither set, then the result depends on if the retain and remove sets are populated.
	 * If retain is populated and the element isn't in it, then return false.
	 * If retain is size 0 and the element isn't in the remove set, return true;
	 * @param element The element to test
	 * @param retain A set of elements to retain
	 * @param remove A set of elements to remove.
	 * @return true if this element should be retained, or false if it should be removed.
	 */
	public boolean retainElement (final T element, final Set<T> retain, final Set<T> remove) {
		if (retain.contains(element)) return true;
		if (remove.contains(element)) return false;
		// element not in retain, but retain has elements.
		if (retain.size()>0) return false;
		// if remove is > 0 OR both sets are empty, return true;
		return true;
	}

}
