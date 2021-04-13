package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.IntervalTagComparator;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;

/**
 * Converts a group of SNP UMI pileups (all the reads that support a single UMI/gene/SNP) into a DigitalAlleleCounts Object, 
 * which tracks the number of UMIs and reads for each expressed allele.
 * 
 * At a single SNP location, there may be more than one gene annotation present.  Reads that overlap multiple genes will be assigned
 * to both genes.  The goal of this iterator is to filter all but one gene from the SNP in a consistent way across cells.  This is done
 * by counting the number of UMIs expressed on each gene across cells, then filtering the SNPUMIBasePileup list to the subset that are 
 * assigned to this gene.
 * @author nemesh
 *
 */
public class DigitalAlleleCountsBestGeneIterator implements DigitalAlleleCountsGeneIteratorI {

	private final GroupingIterator<SNPUMIBasePileup> groupingIter;
	private final SAMSequenceDictionary dict;
	private final int baseQualityThreshold;

	// stash results for a list of cells/genes for a single pileup here.
	private final Queue<DigitalAlleleCounts> stash;

	/**
	 * Construct an iterator that generates DigitalAlleleCount Objects
	 * @param iter An interator across SNPUMIBasePileup objects that contains data about one cell/gene/SNP/UMI for some number of reads
	 * @param baseQualityThreshold A minimum base quality threshold to retain a SNPUMIBasePileup. 
	 */
	public DigitalAlleleCountsBestGeneIterator(SNPUMIBasePileupIterator iter, final int baseQualityThreshold) {
		this.dict = iter.getSNPIntervals().getHeader().getSequenceDictionary();
		// group the data by SNP interval ONLY.

		final Comparator<SNPUMIBasePileup> groupingComparator = new Comparator<SNPUMIBasePileup>() {
			@Override
			public int compare(final SNPUMIBasePileup o1, final SNPUMIBasePileup o2) {
				int ret = IntervalTagComparator.compare(o1.getSNPInterval(), o2.getSNPInterval(), dict);
				return ret;
			}
		};
		this.baseQualityThreshold = baseQualityThreshold;
		this.groupingIter = new GroupingIterator<SNPUMIBasePileup>(iter, groupingComparator);
		this.stash = new LinkedList<>();
	}

	@Override
	public boolean hasNext() {
		if (!stash.isEmpty() || this.groupingIter.hasNext())
			return true;
		return false;
	}

	@Override
	public DigitalAlleleCounts next() {
		if (!stash.isEmpty())
			return stash.poll();
		// no data in stash, populate next SNP.
		populateStash();
		if (!stash.isEmpty())
			return stash.poll();
		// I don't think you can get here unless you don't call "has next".
		return null;
	}

	/**
	 * For a collection of SNPUMIBasePileup objects grouped by SNP to include all cells and genes,
	 * iterator over the collection, find the "best" gene as defined by the gene having the largest number of UMIs, then
	 * filter the collection to only include that gene.  Construct DigitalAlleleCounts objects for each cell, and queue
	 * then so they can be retrieved by calls to next()
	 */
	private void populateStash() {
		// stash is empty, poll the data and build the next SNP interval worth of data.
		// if there's no data left, exit early.
		if (!this.groupingIter.hasNext()) return;
		List<SNPUMIBasePileup> startList = this.groupingIter.next();

		ObjectCounter<String> genePileUpCounter = new ObjectCounter<String>();

		// each SNPUMIBasePileup is the reads for a single UMI.
		for (SNPUMIBasePileup p : startList)
			genePileUpCounter.increment(p.getGene());

		// which gene has the most UMIs across all cells?
		String bestGene = genePileUpCounter.getKeysOrderedByCount(true).get(0);

		Iterator<SNPUMIBasePileup> bestGenePileups = startList.stream().filter(x -> x.getGene().equals(bestGene)).iterator();
		GroupingIterator<SNPUMIBasePileup> groupingIterator = new GroupingIterator<>(bestGenePileups, cellComparator);

		while (groupingIterator.hasNext()) {
			DigitalAlleleCounts dac = DigitalAlleleCountsIterator.getDAC(groupingIterator.next().iterator(), baseQualityThreshold);
			this.stash.add(dac);
		}
	}

	/** 
	 * Group SNPUMIBasePileup objects by cell barcode
	 */
	private static final Comparator<SNPUMIBasePileup> cellComparator = new Comparator<SNPUMIBasePileup>() {
		@Override
		public int compare(final SNPUMIBasePileup o1, final SNPUMIBasePileup o2) {
			int ret = o1.getCell().compareTo(o2.getCell());
			return ret;
		}
	};

	@Override
	public void close() {
		CloserUtil.close(this.groupingIter);
	}

}
