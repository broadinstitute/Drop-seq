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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPUMIBasePileup;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPUMIBasePileupIterator;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SortOrder;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.IntervalTagComparator;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A wrapper around SNPUMIBasePileupIterator to produce SampleGenotypeProbabilities objects.
 *
 * The one "trick" with this iterator is that a SNP pileup can rarely be assigned to multiple genes that we want to avoid double counting.
 * Since we can count on the ordering of snp/gene/cell for records,
 * We group up all pileups for a snp to a list, then filter that list so that each cell is seen once.
 * @author nemesh
 *
 */
public class SampleGenotypeProbabilitiesIterator implements CloseableIterator<SampleGenotypeProbabilities> {

	private final GroupingIterator<SNPUMIBasePileup> groupingIterator;
	private final int editDistance;
	// private final SNPUMIBasePileupIterator groupingIterator;



	/**
	 * Create a new iterator of SampleGenotypeProbabilities
	 * This requires that the SNPUMIBasePileupIterator be sorted in cell order.
	 * This does a few things:
	 * 1) Groups up all pileups for a snp/cell.
	 * 2) If there's more than 1 gene for this snp/cell, pick the gene with the largest number of pileups
	 * 3) Perform edit distance collapse on those pileups, each of which represents a single UMI, and edit distance collapse is
	 * performed on the molecular barode of the pileup.  If you wish to avoid this, set editDistance to 0.
	 *
	 * @param iter A SNPUMIBasePileupIterator that has been generated in sort order: SNPUMIBasePileupIterator.SORT_ORDER.CELL.
	 * @param dict A dictionary file to set
	 * @param editDistance
	 */
	public SampleGenotypeProbabilitiesIterator(final SNPUMIBasePileupIterator iter, final SAMSequenceDictionary dict, final int editDistance, final SortOrder order) {
		// assert that the SNPUMIBasePileupIterator must be sorted in cell before gene order.
		if (iter.getSortOrder()!=order)
			throw new IllegalStateException("When constructing this object, the backing iterator must be sorted in the same order.");
		this.editDistance=editDistance; 

		// group up all cells under a SNP.
		GroupingIterator<SNPUMIBasePileup> groupingIterator=null;
		if (order==SortOrder.SNP_CELL)
			groupingIterator = new GroupingIterator<>(iter,
		            new Comparator<SNPUMIBasePileup>() {
		                @Override
		                public int compare(final SNPUMIBasePileup o1, final SNPUMIBasePileup o2) {
		                	int cmp = IntervalTagComparator.compare(o1.getSNPInterval(), o2.getSNPInterval(), dict);
		                	if (cmp==0)
								cmp = o1.getCell().compareTo(o2.getCell());
		                    return cmp;
		                }
		   });
		else if (order==SortOrder.CELL_SNP)
			groupingIterator = new GroupingIterator<>(iter,
		            new Comparator<SNPUMIBasePileup>() {
		                @Override
		                public int compare(final SNPUMIBasePileup o1, final SNPUMIBasePileup o2) {
		                	int cmp = o1.getCell().compareTo(o2.getCell());
		                	if (cmp==0)
								cmp = IntervalTagComparator.compare(o1.getSNPInterval(), o2.getSNPInterval(), dict);
		                    return cmp;
		                }
		   });
	    this.groupingIterator=groupingIterator;
	}

	public SampleGenotypeProbabilitiesIterator(final SNPUMIBasePileupIterator iter, final SAMSequenceDictionary dict, final int editDistance) {
		this(iter, dict, editDistance, SortOrder.SNP_CELL);
	}

	@Override
	public boolean hasNext() {
		return this.groupingIterator.hasNext();
	}

	@Override
	public SampleGenotypeProbabilities next() {
		return nextPickBestGene();
	}

	/**
	 * Generates a SampleGenotypeProbabilities for a single SNP/cell.
	 * If there are multiple genes, this picks the gene that has the most pileup objects(molecular barcodes)
	 * This reduces a list of pileup objects to a single probabilities object.
	 * This solves the problem where reads could be assigned to multiple genes and the reads were cloned/split from the compound assignment of A/B to two reads with 1 gene each
	 * by selecting the gene with the most support, which gets rid of the "double counting" of reads issue that would otherwise arise, and aggregating results across UMIs of the "best" gene.
	 * @return
	 */
	SampleGenotypeProbabilities nextPickBestGene() {
		// data comes in for a single SNP and cell.
		// want to assign pileups per cell, and pileups for different genes on the same SNP.
		List<SNPUMIBasePileup> startList = groupingIterator.next();
		SNPUMIBasePileup initial = startList.get(0);
		SampleGenotypeProbabilities result = new SampleGenotypeProbabilities(initial.getSNPInterval(), initial.getCell());

		ObjectCounter<String> genePileUpCounter = new ObjectCounter<String>();

		for (SNPUMIBasePileup p: startList)
			genePileUpCounter.increment(p.getGene());
		String bestGene = genePileUpCounter.getKeysOrderedByCount(true).get(0);
		for (SNPUMIBasePileup p: startList)
			if (p.getGene().equals(bestGene))
				result.add(p);
		result.collapseUMIs(this.editDistance);
		return (result);
	}


	/**
	 * Generates a SampleGenotypeProbabilities for a single SNP/cell.
	 * If there are multiple genes, look at the molecular barcodes.
	 * Pick the best pileup (largest number of reads) for each molecular barcode
	 * This reduces a list of pileup objects to a single probabilities object.
	 * @return
	 */
	SampleGenotypeProbabilities nextCombineGenes() {
		// data comes in for a single SNP and cell.
		// want to assign pileups per cell, and pileups for different genes on the same SNP.
		List<SNPUMIBasePileup> startList = groupingIterator.next();

		Map<String, SNPUMIBasePileup> bestPileUpForMolBCMap = new HashMap<String, SNPUMIBasePileup>();

		for (SNPUMIBasePileup currentPileUp: startList) {
			String molBC = currentPileUp.getMolecularBarcode();
			SNPUMIBasePileup bestP = bestPileUpForMolBCMap.get(molBC);
			// if there's no pileup for this molecular barcode, assign this as best.
			if (bestP==null)
				bestPileUpForMolBCMap.put(molBC, currentPileUp);
			else { // otherwise, test the best against the current.
				int numReadsBest = bestP.getNumBases();
				int numReadsCurrent = currentPileUp.getNumBases();
				if (numReadsCurrent>numReadsBest)
					bestPileUpForMolBCMap.put(molBC, currentPileUp);
			}
		}

		// set up result.
		SNPUMIBasePileup initial = startList.get(0);
		SampleGenotypeProbabilities result = new SampleGenotypeProbabilities(initial.getSNPInterval(), initial.getCell());

		// add pileups to result.
		for (SNPUMIBasePileup p: bestPileUpForMolBCMap.values())
			result.add(p);
		result.collapseUMIs(this.editDistance);
		return (result);
	}


	@Override
	public void remove() {
		this.groupingIterator.remove();
	}

	@Override
	public void close() {
        CloserUtil.close(this.groupingIterator);
	}




}
