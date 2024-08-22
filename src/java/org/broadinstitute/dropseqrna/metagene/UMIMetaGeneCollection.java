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

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.dropseqrna.metagene.MetaGene;
import org.broadinstitute.dropseqrna.metagene.ReadGroupResult.MetaGeneTypeEnum;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.ReadNameComparator;
import org.broadinstitute.dropseqrna.utils.readiterators.MapQualityPredicate;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * For all the reads on a cell/UMI, discover the meta-genes.
 *
 * There are a number of cases that are handled here:
 * 1) Reads only map to a single gene uniquely - there's no meta gene or ambiguous info
 * 2) Reads map to a single gene uniquely, and to more than one gene (for reads that map non-uniquely).  In this case
 * the reads are all assigned to the unique gene.
 * 3) Meta gene case 1 (unambiguous): reads map non-uniquely to 2 genes, but not to either gene uniquely.  These reads will be considered a meta-gene.
 * 4) Meta gene case 2 (ambiguous): reads map non-uniquely to 2 genes, but also to both genes uniquely (as opposed to case 2 above).
 * This case may or may not happen, and these reads are either errors, or something incredibly interesting is happening.
 * These cases will be tracked.
 *
 * A gene can only be in one category: unique, unambiguous meta-gene, or ambiguous meta-gene.
 *
 * This class will also take a list of reads, and given a map of approved meta genes (one gene -> set of all related genes in meta gene), will tag meta gene annotations onto the reads appropriately.
 * In this case, reads that map uniquely are left alone.  Reads that map to a gene in the list and the other reads map to genes in the set are tagged as meta genes.  Reads that map to a gene in the set
 * and another gene not in the set are left as a unique gene.
 *
 *
 *
 * @author nemesh
 *
 */



//Set<Set<String>> instead of Set<List<String>>
//create a class that is comparable and set and is immutable.
//or use ImmutableSortedSet


public class UMIMetaGeneCollection {

	private final int mapQualityUnique;
	private final String geneTag;
	private final String molBCTag;
	private final String cellBCTag;
	private final String cellBarcode;
	private final String molecularBarcode;
	private final List<SAMRecord> informativeReads;
	private Set<String> uniqueGenes;
	private Set<MetaGene> ambiguousMetaGenes;
	private Set<MetaGene> unambiguousMetaGenes;

	private static final Log log = Log.getInstance(UMIMetaGeneCollection.class);
	// last time you asked for unique genes, etc, did you use edit distance on the reads?
	// if this agrees, then you can use cached results.
	// private boolean curentUseEditDistance=false;

	public UMIMetaGeneCollection (final String geneTag, final String molBCTag, final String cellBCTag, final Collection<SAMRecord> informativeReads, final int mapQualityUnique) {
		this.mapQualityUnique=mapQualityUnique;
		this.geneTag=geneTag;
		this.molBCTag=molBCTag;
		this.cellBCTag=cellBCTag;
		List<SAMRecord> infRecs = new ArrayList<>(informativeReads);
		Collections.sort(infRecs, new SAMRecordQueryNameComparator());
		this.informativeReads=infRecs;		
		this.uniqueGenes=null;
		this.ambiguousMetaGenes=null;
		this.unambiguousMetaGenes=null;
		if (informativeReads.size()>0) {
			this.cellBarcode=infRecs.get(0).getStringAttribute(cellBCTag);
			this.molecularBarcode=infRecs.get(0).getStringAttribute(molBCTag);
		} else {
			this.cellBarcode=null;
			this.molecularBarcode=null;
		}
	}

	/**
	 * For the informative reads for this UMI/Cell, group the informative reads by queryname, and decide which groups are unique/metagenes/etc.
	 * @return
	 */
	public Collection<ReadGroupResult> getResultsPerRead () {
		List<ReadGroupResult> result = new ArrayList<>();
		// Set<String> uniqueGenes = this.getUniquelyMappedGenes();
		GroupingIterator<SAMRecord> gIter = new GroupingIterator<>(informativeReads.iterator(), new ReadNameComparator());
		for (List<SAMRecord> recs: gIter) {
			ReadGroupResult rgr = getMetaGeneType(recs);
			result.add(rgr);
		}
		
		try {
			gIter.close();
		} catch (IOException e) {
			// this should never happen.
			e.printStackTrace();
		}				
		return (result);
	}



	/**
	 * Unambiguous genes are those where reads for this UMI only map to a single gene uniquely.
	 * If a read maps to a gene uniquely
	 * @return
	 */
	public Set<MetaGene> getUnambiguousMetaGenes() {
		if (this.unambiguousMetaGenes!=null)
			return this.unambiguousMetaGenes;
		populateGeneSets();
		return this.unambiguousMetaGenes;
	}

	public Set<MetaGene> getAmbiguousMetaGenes() {
		if (this.ambiguousMetaGenes!=null)
			return this.ambiguousMetaGenes;
		populateGeneSets();
		return this.ambiguousMetaGenes;
	}

	/**
	 * Get the set of genes where reads for this UMI only map to a single gene.
	 * These are the uniquely mapped genes.
	 * @return
	 */
	public Set<String> getUniquelyMappedGenes () {
		// GroupingIterator<SAMRecord> gIter = new GroupingIterator<>(reads.iterator(), new ReadNameComparator());
		if (uniqueGenes!=null) return this.uniqueGenes;
		this.uniqueGenes = new HashSet<>();

		MapQualityPredicate p = new MapQualityPredicate(this.mapQualityUnique, false);
		for (SAMRecord r: informativeReads) {
			if (p.test(r)) {
				String gene = r.getStringAttribute(this.geneTag);
				this.uniqueGenes.add(gene);
			}
		}

		return this.uniqueGenes;
	}

	/**
	 * For the informative reads, group by read name, get the best read(s) for each read name.
	 * If there's only one, it's a "unique" mapping.
	 * This promotes reads with low map quality but better edit distance to look like unique mappings.
	 *
	 * @return
	 */
	private Set<String> getUniquelyMappedGenesWithEDComparison() {
		this.uniqueGenes = new HashSet<>();
		GroupingIterator<SAMRecord> gIter = new GroupingIterator<>(informativeReads.iterator(), new ReadNameComparator());

		while (gIter.hasNext()) {
			List<SAMRecord> recs = gIter.next();
			List<SAMRecord> recsUnique = new ArrayList<>(getBestReadsByEditDistance(recs));
			if (recsUnique.size()==1) {
				SAMRecord r = recsUnique.get(0);
				String gene = r.getStringAttribute(this.geneTag);
				this.uniqueGenes.add(gene);
			}
		}
		CloserUtil.close(gIter);
		return this.uniqueGenes;
	}


	/**
	 *
	 */
	void populateGeneSets () {
		
		Set<String> uniqueGenes = this.getUniquelyMappedGenes();
		Set<MetaGene> ambiguous = new HashSet<>();
		Set<MetaGene> unambiguous = new HashSet<>();

		GroupingIterator<SAMRecord> gIter = new GroupingIterator<>(informativeReads.iterator(), new ReadNameComparator());

		// process data by read name.  One read name with 1 or more reads is a ReadGroupResult
		for (List<SAMRecord> list : gIter) {
			ReadGroupResult r = getMetaGeneType (list);
			MetaGene mg = new MetaGene(r.getGenes());
			if (r.getType().equals(ReadGroupResult.MetaGeneTypeEnum.UNAMBIGIOUS))
				unambiguous.add(mg);
			if (r.getType().equals(ReadGroupResult.MetaGeneTypeEnum.AMBIGUOUS))
				ambiguous.add(mg);
		}
		// these are mutually exclusive.
		this.unambiguousMetaGenes=unambiguous;
		this.ambiguousMetaGenes=ambiguous;
		// pretty well validated, no need to waste time on the call.
		// validate();
		
		try {
			gIter.close();
		} catch (IOException e) {
			// This should never, ever throw an error.
			e.printStackTrace();
		}		
	}

	/**
	 * Given some reads for this UMI have already been mapped uniquely to genes and share a read name, decide the state of these mappings.
	 * 
	 * Basic algorithm pseudocode:
	 * 1) Get the unique gene names for this set of reads [this can be less than the number of reads in the list if reads are unmapped]
	 * 2) If there's a single read, classify as unmapped [map quality=0], low MQ [low map quality], high map quality [unique] or other [missing classification status]
	 * 3) If there's more than 1 read, and that read maps to multiple genes:
	 * 		If 1 of the genes is uniquely mapped elsewhere, then this read supports the uniquely mapped one.
	 * 		If none of the genes is uniquely mapped elsewhere, then this read is a meta gene.
	 * 		If these reads map to multiple genes that map uniquely elsewhere, then this read is ambiguous.
	 * @param list A list of SAMRecords to look at within the context of a Cell/UMI.
	 * @return
	 */
	public ReadGroupResult getMetaGeneType (final List<SAMRecord> list) {
		String readName = list.get(0).getReadName();
		Set<String> readGenes = new HashSet<>();

		for (SAMRecord r: list) {
			if (!r.getReadName().equals(readName)) throw new IllegalStateException("Read names don't match for reads that should all share the same name");
			String gene = r.getStringAttribute(this.geneTag);
			if (gene!=null) readGenes.add(gene);
		}

		List<String> orderedReadGenes = new ArrayList<>(readGenes);

		// if there's a single read, the read is unmapped or unique, or I'm dumb and missed something (OTHER)
		if (list.size()<2) {
			int readMQ = list.get(0).getMappingQuality();
			if (readMQ==0) return new ReadGroupResult(orderedReadGenes, MetaGeneTypeEnum.UNMAPPED, list);
			if (readMQ<this.mapQualityUnique) return new ReadGroupResult(orderedReadGenes, MetaGeneTypeEnum.LOW_MQ, list);
			if (readMQ>=this.mapQualityUnique) return new ReadGroupResult(orderedReadGenes, MetaGeneTypeEnum.UNIQUE, list);;
			return new ReadGroupResult(orderedReadGenes, MetaGeneTypeEnum.OTHER, list);
		}


		Set<String> intersect = new HashSet<>(readGenes);
		// set intersect - if the possible locations all map to 1 gene that is mapped uniquely, this is not a meta gene.
		// if all possible locations map to more than 1 gene that is also uniquely mapped, these genes are ambiguous.
		intersect.retainAll(uniqueGenes);

		// the mapping to the unique gene overrides the fact that the read maps to >1 gene.
		if (readGenes.size()>1 & intersect.size()==1)
			return new ReadGroupResult(intersect, MetaGeneTypeEnum.UNIQUE, list);

		// unambiguous meta gene - no unique mappings of the genes.
		if (readGenes.size()>1 & intersect.size()==0)
			return new ReadGroupResult(readGenes, MetaGeneTypeEnum.UNAMBIGIOUS, list);

		// ambiguous meta gene - both genes have unique mappings, and this read group shares both mappings.
		if (readGenes.size()>1 & intersect.size()>1)
			return new ReadGroupResult(readGenes, MetaGeneTypeEnum.AMBIGUOUS, list);

		// covering the case that I didn't cover a case...
		return new ReadGroupResult(readGenes, MetaGeneTypeEnum.OTHER, list);
	}

	/**
	 * Generate the result for this (possibly) multi-read given the set of genes that uniquely map.
	 * @param uniquelyMappedGenes
	 * @param list
	 * @param useEditDistance
	 * @return
	 */ 
	/*
	//TODO: this is terrible cut and paste and should be removed.
	public ReadGroupResult getMetaGeneType (final Set<String> uniquelyMappedGenes, final List<SAMRecord> list, final Boolean useEditDistance) {
		if (useEditDistance==null || useEditDistance==false) return getMetaGeneType(list);

		String readName = list.get(0).getReadName();
		Set<String> readGenes = new HashSet<>();

		// pick out the read(s) with the best edit distances.
		Collection<SAMRecord> filtered = getBestReadsByEditDistance(list);

		// what genes do I have in the filtered data set?
		for (SAMRecord r: filtered) {
			if (!r.getReadName().equals(readName)) throw new IllegalStateException("Read names don't match for reads that should all share the same name");
			String gene = r.getStringAttribute(this.geneTag);
			if (gene!=null) readGenes.add(gene);
		}

		List<String> orderedReadGenes = new ArrayList<>(readGenes);

		// if there's a single read after filtering, the read is unmapped or unique.
		if (list.size()<2 && filtered.size()<2) {
			if (list.get(0).getReadUnmappedFlag()) return new ReadGroupResult(orderedReadGenes, MetaGeneTypeEnum.UNMAPPED, list);
			return new ReadGroupResult(orderedReadGenes, MetaGeneTypeEnum.UNIQUE, list);
		}

		// since we've filtered reads by edit distance, there can be 1 or more reads, but they are the best edit distance.
		Set<String> intersect = new HashSet<>(readGenes);

		// set intersect - if the possible locations all map to 1 gene that is mapped uniquely, this is not a meta gene.
		// if all possible locations map to more than 1 gene that is also uniquely mapped, these genes are ambiguous.
		intersect.retainAll(uniqueGenes);

		// unambiguous meta gene - no unique mappings of the genes.
		if (readGenes.size()>1 & intersect.size()==0)
			return new ReadGroupResult(readGenes, MetaGeneTypeEnum.UNAMBIGIOUS, list);

		// ambiguous meta gene - both genes have unique mappings, and this read group shares both mappings.
		if (readGenes.size()>1 & intersect.size()>1)
			return new ReadGroupResult(readGenes, MetaGeneTypeEnum.AMBIGUOUS, list);

		// covering the case that I didn't cover a case...
		return new ReadGroupResult(readGenes, MetaGeneTypeEnum.OTHER, list);
	}
	*/
	
	/**
	 * Get the read(s) that have the lowest edit distance.  As long as there is at least one read in the input, one read will be returned.
	 * @param recs
	 * @return
	 */
	public static Collection<SAMRecord> getBestReadsByEditDistance (final Collection<SAMRecord> recs) {
		if (recs==null) return null;
		if (recs.size()<2) return recs;
		// find the minimum edit distance of the reads.
		int bestED = recs.stream().filter(x -> x.getIntegerAttribute("NM")!=null).mapToInt(r -> r.getIntegerAttribute("NM")).min().getAsInt();
		// find the reads at that minimum edit distance.
		List<SAMRecord> finalRecs = recs.stream().filter(x -> x.getIntegerAttribute("NM")!=null).filter(r -> r.getIntegerAttribute("NM") <= bestED).collect(Collectors.toList());
		return finalRecs;
	}

	/**
	 * A gene can only be classified one way: Ambiguous meta gene or Unambiguous meta gene
	 * A UMI's unique genes should be a non-overlapping set with it's unambiguous meta genes.
	 * @return
	 */
	public boolean validate() {
		// do the ambiguous meta genes and unambigous meta genes overlap?
		Set<MetaGene> test1 = new HashSet<>(this.unambiguousMetaGenes);
		test1.retainAll(this.ambiguousMetaGenes);
		if (test1.size()>0) {
			log.error("UMI contains genes that are both ambiguous and unambiguous meta genes "+ this.toString());
			return false;
		}
		// are any unique genes assigned to unambigous meta genes for this UMI?
		Set<MetaGene> test2 = new HashSet<>(this.unambiguousMetaGenes);
		test2.retainAll(this.uniqueGenes);
		if (test2.size()>0) {
			log.error("UMI contains genes that are both unambiguous and unique " + this.toString());
			return false;
		}
		return true;
	}

	public String getMolBCTag() {
		return molBCTag;
	}
 
	public String getCellBCTag() {
		return cellBCTag;
	}

	public String getCellBarcode() {
		return cellBarcode;
	}

	public String getMolecularBarcode() {
		return molecularBarcode;
	}

	public int getNumInformativeReads () {
		return this.informativeReads.size();
	}

	@Override
	public String toString () {
		String r = "Cell Barcode [" +this.cellBarcode+"] molBC ["+ this.molecularBarcode+"] num reads [" + this.informativeReads.size() +"] unique genes "+ this.getUniquelyMappedGenes() +" meta genes " + this.getUnambiguousMetaGenes().toString()+" ambiguous metagenes " + this.getAmbiguousMetaGenes();
		return r;
	}



}
