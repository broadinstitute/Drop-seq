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
package org.broadinstitute.dropseqrna.metagene;

import com.google.common.collect.ImmutableSortedSet;
import org.apache.commons.lang3.StringUtils;

import java.util.Collection;
import java.util.Iterator;

/**
 * Holds a set of gene names that form a meta gene.
 * This object is immutable, ordered, and comparable.
 * @author nemesh
 *
 */
public class MetaGene implements Comparable<MetaGene> {

	private final ImmutableSortedSet<String> geneSet;

	public MetaGene (final Collection<String> geneNames) {
		geneSet = ImmutableSortedSet.copyOf(geneNames);
	}

	public ImmutableSortedSet<String> getGeneNames () {
		return this.geneSet;
	}

	/**
	 * Sort all the gene names, then concatenate to a single string separated by the sep character.
	 * @param sep
	 * @return
	 */
	public String getMetaGeneName (final char sep) {
		return StringUtils.join(this.geneSet, sep);
	}
	
	/**
	 * Are all the genes in the other object contained in this metagene object
	 * For example, if other contains A, and this object contains A,B, this would return true.
	 * @param other A metagene to compare against this for a partial match
	 * @return True if this metagene contains all genes from the other object.
	 */
	public boolean partialMatch (MetaGene other) {
		for (String otherGene: other.geneSet) {
			if (!geneSet.contains(otherGene)) return false;
		}
		return true;
	}
	
	@Override
	public int compareTo(final MetaGene o) {
		int result =0;
		Iterator<String> iter1= geneSet.iterator();
		Iterator<String> iter2= o.geneSet.iterator();
		while (iter1.hasNext() && iter2.hasNext()) {
			String s1 = iter1.next();
			String s2 = iter2.next();
			result = s1.compareTo(s2);
			if (result!=0) return result;
		}
		// results equal so far.  smaller sets come first.
		result= geneSet.size() - o.geneSet.size();
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
		MetaGene other = (MetaGene) obj;
		if (geneSet == null) {
            return other.geneSet == null;
		} else return geneSet.equals(other.geneSet);
    }

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((geneSet == null) ? 0 : geneSet.hashCode());
		return result;
	}

	@Override
	public String toString() {
		return getMetaGeneName(':');
	}

}
