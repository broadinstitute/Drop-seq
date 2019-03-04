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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.PeekableIterator;

/**
 * Constructs an iterator to produce MultiCellDigitalAlleleCounts objects from a bunch of DigitalAlleleCounts objects. Iteratorception-ception!
 * A MultiCellDigitalAlleleCounts object is just the aggregation of all DigitalAlleleCounts objects across all data for a snp/gene.
 * @author nemesh
 *
 */

public class MultiCellDigitalAlleleCountsIterator implements CloseableIterator<MultiCellDigitalAlleleCounts> {

	private final PeekableIterator<DigitalAlleleCounts> iter;
	
	public MultiCellDigitalAlleleCountsIterator (DigitalAlleleCountsIterator iter) {
		this.iter=new PeekableIterator<DigitalAlleleCounts>(iter);
	}
	
	
	@Override
	public MultiCellDigitalAlleleCounts next() {
		if (!iter.hasNext()) return (null);
		
		// get the pileup object and make the initial DAC with the first pileup.
		DigitalAlleleCounts dac = iter.next();
		
		MultiCellDigitalAlleleCounts multiDAC = new MultiCellDigitalAlleleCounts(dac.getSnpInterval(), dac.getGene());
		String currentGene = dac.getGene();
		Interval currentSNP = dac.getSnpInterval();
		multiDAC.add(dac);
		
		// get SNPUMIBasePileup objects until the gene, cell, or snpInterval changes.
		while (this.iter.hasNext()) {
			dac=iter.peek(); // check the next object to see if it should be added.
			if (dac.getGene().equals(currentGene) && dac.getSnpInterval().equals(currentSNP)) {
				dac=iter.next(); // convert from the peeked object to the real object to advance the iterator.
				multiDAC.add(dac);
			} else {
				 break; // the next peeked object was a a pileup for a different DAC, quit the loop.
			}
		}		
		return (multiDAC);
				
		
	}

	@Override
	public void remove() {
		this.iter.remove();
	}

	@Override
	public void close() {
		CloserUtil.close(this.iter);
	}

	@Override
	public boolean hasNext() {
		return this.iter.hasNext();
	}
	
	
	
}
