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

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.IntervalTagComparator;

import java.util.Comparator;
import java.util.Iterator;

/**
 * Constructs an iterator to produce DigitalAlleleCounts objects from a SNPUMIBasePileupIterator. Iteratorception!
 * A DigitalAlleleCounts object is just the sum of all SNPUMIBasePileup objects across all molecular barcodes for a snp/cell/gene.
 * @author nemesh
 *
 */
public class DigitalAlleleCountsIterator implements CloseableIterator<DigitalAlleleCounts>{

	private final GroupingIterator<SNPUMIBasePileup> groupingIter;

	private final int baseQualityThreshold;
	private final SAMSequenceDictionary dict;

	public DigitalAlleleCountsIterator (final SNPUMIBasePileupIterator iter, final int baseQualityThreshold) {

		this.baseQualityThreshold=baseQualityThreshold;
		this.dict=iter.getSNPIntervals().getHeader().getSequenceDictionary();

        final Comparator<SNPUMIBasePileup> groupingComparator = new Comparator<SNPUMIBasePileup>() {
            @Override
            public int compare(final SNPUMIBasePileup o1, final SNPUMIBasePileup o2) {
            	int ret = IntervalTagComparator.compare(o1.getSNPInterval(), o2.getSNPInterval(), dict);
                //int ret = o1.getSnpID().compareTo(o2.getSnpID());
                if (ret == 0) {
                    ret = o1.getGene().compareTo(o2.getGene());
                    if (ret == 0) {
                        ret = o1.getCell().compareTo(o2.getCell());
                    }
                }
                return ret;
            }
        };
		this.groupingIter = new GroupingIterator<SNPUMIBasePileup>(iter, groupingComparator);
	}

	/**
	 * Return the next DigitalAlleleCounts object.
	 *
	 * @return The next DAC object, or null if the iterator is exhausted.
	 */
	@Override
	public DigitalAlleleCounts next() {
		if (!groupingIter.hasNext()) return (null);

		// get the pileup object and make the initial DAC with the first pileup.
		final Iterator<SNPUMIBasePileup> pileupIter = groupingIter.next().iterator();
        DigitalAlleleCounts dac = getDAC(pileupIter, this.baseQualityThreshold); 		
		return dac;
	}
	
	public static DigitalAlleleCounts getDAC (Iterator<SNPUMIBasePileup> pileupIter, int baseQualityThreshold) {
		SNPUMIBasePileup p = pileupIter.next();

		final DigitalAlleleCounts dac = new DigitalAlleleCounts(p.getSNPInterval(), p.getGene(), p.getCell(), baseQualityThreshold);
		dac.addPileup(p);

		// get SNPUMIBasePileup objects until the gene, cell, or snpInterval changes.
		while (pileupIter.hasNext()) {
				dac.addPileup(pileupIter.next());
		}
		return dac;
	}

	@Override
	public void remove() {
		this.groupingIter.remove();
	}

	@Override
	public void close() {
		CloserUtil.close(this.groupingIter);
	}

	@Override
	public boolean hasNext() {
		return this.groupingIter.hasNext();
	}
}
