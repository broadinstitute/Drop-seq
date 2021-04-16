package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;

import org.apache.commons.lang.builder.CompareToBuilder;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.IntervalTagComparator;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;



import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;

/**
 * Converts a group of SNP UMI pileups (all the reads that support a single UMI/gene/SNP) into a DigitalAlleleCounts Object, 
 * which tracks the number of UMIs and reads for each expressed allele.
 * 
 * At a single SNP location, there may be more than one gene annotation present.  Reads that overlap multiple genes will be assigned
 * to both genes.  The goal of this iterator is to filter all but one gene from the SNP in a consistent way across cells.  This is done
 * by counting the number of UMIs expressed on each gene across cells, then filtering the SNPUMIBasePileup list to the subset that are 
 * assigned to this gene.
 * 
 * This is heavily modeled after CountChangingIteratorWrapper's method of having a queue and processing it.
 * @author nemesh
 *
 */
public class DigitalAlleleCountsBestGeneIterator implements DigitalAlleleCountsGeneIteratorI {

	private final GroupingIterator<SNPUMIBasePileup> groupingIter;
	private final SAMSequenceDictionary dict;
	private final int baseQualityThreshold;
	private final Map<Interval, String> refAllele;
	private final Map<Interval, String> altAllele;
	
	// stash results for a list of cells/genes for a single pileup here.
	private final Queue<DigitalAlleleCounts> stash;

	/**
	 * Construct an iterator that generates DigitalAlleleCount Objects
	 * @param iter An interator across SNPUMIBasePileup objects that contains data about one cell/gene/SNP/UMI for some number of reads
	 * @param baseQualityThreshold A minimum base quality threshold to retain a SNPUMIBasePileup. 
	 */
	public DigitalAlleleCountsBestGeneIterator(SNPUMIBasePileupIterator iter, final int baseQualityThreshold, Map<Interval, String> refAllele, Map<Interval, String> altAllele) {
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
		this.refAllele=refAllele;
		this.altAllele=altAllele;
	}

	@Override
	public boolean hasNext() {
		populateStash();
        return !stash.isEmpty();        	
	}

	@Override
	public DigitalAlleleCounts next() {
		populateStash();
        return stash.remove();        		
	}
	
	private void populateStash() {
		while (stash.isEmpty() && this.groupingIter.hasNext()) {
			processRecordGroup(this.groupingIter.next());
        }				
	}
	
	/**
	 * For a collection of SNPUMIBasePileup objects grouped by SNP to include all cells and genes,
	 * iterator over the collection, find the "best" gene as defined by the gene having the largest number of UMIs, then
	 * filter the collection to only include that gene.  Construct DigitalAlleleCounts objects for each cell, and queue
	 * then so they can be retrieved by calls to next()
	 */	
	private void processRecordGroup (List<SNPUMIBasePileup> startList) {
		// stash is empty, poll the data and build the next SNP interval worth of data.				
				
		if (startList.get(0).getPosition()==154823512)
			System.out.println("STOP");
								
		List<GeneScore> scores = getGeneScoreList  (startList);
		String bestGene=scores.get(0).name;
				
		Iterator<SNPUMIBasePileup> bestGenePileups = startList.stream().filter(x -> x.getGene().equals(bestGene)).iterator();
		GroupingIterator<SNPUMIBasePileup> groupingIterator = new GroupingIterator<>(bestGenePileups, cellComparator);

		while (groupingIterator.hasNext()) {
			DigitalAlleleCounts dac = DigitalAlleleCountsIterator.getDAC(groupingIterator.next().iterator(), baseQualityThreshold, this.refAllele, this.altAllele);
			this.stash.add(dac);
		}
	}
	
	private List<GeneScore> getGeneScoreList  (List<SNPUMIBasePileup> pileups) {
		Map <String, GeneScore> map = new HashMap<>();
		for (SNPUMIBasePileup p: pileups) {
			GeneScore s  = map.get(p.getGene());
			if (s==null) {
				s = new GeneScore(p.getGene());
				map.put(p.getGene(), s);
			}
			s.add(p);
		}
		
		List<GeneScore> scores = new ArrayList<> (map.values());		
		// 0th element best element.
		Collections.sort(scores, Collections.reverseOrder());				
		return (scores);		
		
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
	
	/**
	 * Need to track multiple features of a gene across many SNPUMIBasePileup objects to generate the "best" pileup.
	 * @author nemesh
	 *
	 */
	private class GeneScore implements Comparable <GeneScore>{
		private final String name;
		private int umiCount;
		private Mean averageBaseQuality;
		private double umiPurity;
		
		public GeneScore (String name) {
			this.name = name;
			averageBaseQuality = new Mean();
		}
		
		public void add (SNPUMIBasePileup p) {
			if (!p.getGene().equals(this.name)) 
				throw new IllegalArgumentException("Adding the wrong gene to this struct");
			umiCount++;
			p.getQualities().stream().forEach(x -> this.averageBaseQuality.increment(x));
			this.umiPurity=p.getUMIPurity();
		}
		
		public int getUmiCount () {
			return this.umiCount;
		}
		
		public String getGene () {
			return this.name;
		}
		
		public int getNumReads () {
			return (int) this.averageBaseQuality.getN();
		}
		
		public Double getAverageGenotypeQuality () {
			return this.averageBaseQuality.getResult();			
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;			
			result = prime * result + ((getAverageGenotypeQuality() == null) ? 0 : getAverageGenotypeQuality().hashCode());
			result = prime * result + ((name == null) ? 0 : name.hashCode());
			result = prime * result + umiCount;
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			GeneScore other = (GeneScore) obj;
			if (getAverageGenotypeQuality() == null) {
				if (other.getAverageGenotypeQuality() != null)
					return false;
			} else if (!getAverageGenotypeQuality().equals(other.getAverageGenotypeQuality()))
				return false;
			if (name == null) {
				if (other.name != null)
					return false;
			} else if (!name.equals(other.name))
				return false;
			if (umiCount != other.umiCount)
				return false;
			return true;
		}

		@Override
		public int compareTo(GeneScore o) {
			return new CompareToBuilder()	
						.append(this.umiPurity, o.umiPurity)
				       .append(this.umiCount, o.umiCount)
				       .append(this.getAverageGenotypeQuality(), this.getAverageGenotypeQuality())
				       .append(this.getNumReads(), o.getNumReads())				       
				       .toComparison();				       			
		}

		@Override
		public String toString () {
			StringBuilder b = new StringBuilder();
			b.append("Gene ["+this.name+"] ");
			b.append("UMIs [" + this.umiCount+"] ");
			b.append("Average BQ [" + this.getAverageGenotypeQuality() +"] ");
			b.append("Num Reads [" + this.getNumReads() +"] ");
			b.append("UMI Purity ["+ this.umiPurity+"]");
			return b.toString();
		}
		
 	}
	
	

}
