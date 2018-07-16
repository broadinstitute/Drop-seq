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

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.junit.Assert;
import org.testng.annotations.Test;

public class RetainRemoveListTest {

	@Test
	public void testRemoveOnly() {
		String [] elements = {"1", "2", "3", "5"};
		String [] remove = {"1", "4"};
		String [] expected = {"2", "3", "5"};
		RetainRemoveList<String> rrl = new RetainRemoveList<>();
		List<String> retainedElements = rrl.getElementsToRetain(Arrays.asList(elements), getArrayAsSet(remove), new HashSet<String>());
		String[] re = new String[retainedElements.size()];
		re = retainedElements.toArray(re);
		Assert.assertArrayEquals(expected, re);

	}

	@Test
	public void testRetainOnly() {
		String [] elements = {"1", "2", "3", "5"};
		String [] retain = {"1", "4", "5"};
		String [] expected = {"1", "5"};
		RetainRemoveList<String> rrl = new RetainRemoveList<>();
		List<String> retainedElements = rrl.getElementsToRetain(Arrays.asList(elements), new HashSet<String>(), getArrayAsSet(retain));
		String[] re = new String[retainedElements.size()];
		re = retainedElements.toArray(re);
		Assert.assertArrayEquals(expected, re);
	}

	@Test
	public void testRetainRemove() {
		String [] elements = {"1", "2", "3", "5"};
		String [] retain = {"1", "4"};
		String [] remove = {"1", "5"};
		// retain overrides remove, so we expect to not remove 1.
		String [] expected = {"1"};
		RetainRemoveList<String> rrl = new RetainRemoveList<>();
		List<String> retainedElements = rrl.getElementsToRetain(Arrays.asList(elements), getArrayAsSet(remove), getArrayAsSet(retain));
		String[] re = new String[retainedElements.size()];
		re = retainedElements.toArray(re);
		Assert.assertArrayEquals(expected, re);

	}

	public Set<String> getArrayAsSet (final String [] elements) {
		return new HashSet<> (Arrays.asList(elements));
	}


}
