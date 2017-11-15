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
package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.Iterator;

import org.broadinstitute.dropseqrna.barnyard.Utils;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;

import htsjdk.samtools.SAMRecord;

/**
 * Filters reads where the read strand and gene strand disagree.
 * @author nemesh
 *
 */
public class GeneStrandFilteringIterator extends FilteredIterator<SAMRecord>{

	private final String strandTag;

	protected GeneStrandFilteringIterator(final Iterator<SAMRecord> underlyingIterator, final String strandTag) {
		super(underlyingIterator);
		this.strandTag=strandTag;
	}

	@Override
	/**
	 * For a given record, if the gene strand tag does not match the read strand tag, reject the record.
	 * If there are multiple gene strand tags, they must all agree with the read strand tag.
	 * @param rec
	 * @return
	 */
	public boolean filterOut(final SAMRecord rec) {
		/*
		if (rec.getReadName().equals("HN7TNBGXX:4:11609:2753:3252"))
			System.out.println("STOP");
		*/
		String geneStrand = rec.getStringAttribute(strandTag);
		if (geneStrand==null) return false;

		String [] strands = geneStrand.split(",");
		String readStrandString = Utils.strandToString(!rec.getReadNegativeStrandFlag());
		for (String s: strands)
			if (!s.equals(readStrandString)) return true;
		return false;
	}


}
