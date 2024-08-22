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

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;

import org.apache.commons.lang3.StringUtils;

import htsjdk.samtools.SAMRecord;
import org.broadinstitute.dropseqrna.metagene.MetaGene;

/**
 * For one or more mappings of a read, what genes is the read mapped to, and what type of meta gene (if any) is this?
 * @author nemesh
 *
 */
public class ReadGroupResult {

	private final List<String> genes;
	private MetaGene metaGene;
	
	private final MetaGeneTypeEnum type;
	private final List<SAMRecord> reads;
	
	public ReadGroupResult(final Collection<String> genes, final MetaGeneTypeEnum type, final List<SAMRecord> reads) {
		List<String> g = new ArrayList<>(genes);
		Collections.sort(g);
		this.genes=g;
		this.type=type;
		this.reads=reads;
		this.metaGene=new MetaGene(this.genes);
	}
	
	public List<String> getGenes() {
		return genes;
	}
	
	public MetaGene getMetaGene () {
		return this.metaGene;
	}
	
	/**
	 * For a given read tag, get the unique, sorted set of values.
	 * @return
	 */
	/*
	public Set<String> getTagValue (String tag) {
		Set<String> s =reads.stream().map(x -> x.getStringAttribute(tag)).collect(Collectors.toSet());
		SortedSet<String> values = new TreeSet<>(s) ;
		return values;		
	}
	*/
	
	public MetaGeneTypeEnum getType() {
		return type;
	}

	public enum MetaGeneTypeEnum {
		UNIQUE,UNAMBIGIOUS,AMBIGUOUS,UNMAPPED,LOW_MQ,SPLIT,OTHER
	}

	public String getMetaGeneName (final char sep) {
		return StringUtils.join(genes, sep);
	}

	@Override
	public String toString () {
		return "Genes [" + StringUtils.join(genes, ',') +"] Type ["+ type +"]";
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((genes == null) ? 0 : genes.hashCode());
		result = prime * result + ((type == null) ? 0 : type.hashCode());
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
		ReadGroupResult other = (ReadGroupResult) obj;
		if (genes == null) {
			if (other.genes != null)
				return false;
		} else if (!genes.equals(other.genes))
			return false;
		if (type != other.type)
			return false;
		return true;
	}

	public List<SAMRecord> getReads() {
		return reads;
	}

}
