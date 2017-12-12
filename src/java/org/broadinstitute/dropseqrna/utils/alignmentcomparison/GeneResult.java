/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.dropseqrna.utils.alignmentcomparison;

import java.util.Collection;

import org.broadinstitute.dropseqrna.utils.ObjectCounter;

import com.sun.jdi.request.InvalidRequestStateException;

import htsjdk.samtools.util.StringUtil;

/**
 * For each gene, how many of the reads uniquely mapped to the first gene mapped to the same gene in the second data set
 * Where else did reads align to uniquely?
 * Where else did reads align non-uniquely?
 * @author nemesh
 *
 */
public class GeneResult {

	private String originalGene;
	private String originalContig;
	private int countOriginalReads;
	private int countSameMapping;

	// count of reads that map to a different area uniquely.
	// can I break this up into reads that map to genes and reads that don't map at all?
	private int countDifferentUnique;
	private int countNoGene;

	// the original gene mapping is present but read also maps to some other non-gene region, making the map quality of the gene read low so it wouldn't be in the DGE.
	// this read only contains 1 gene.
	private int countMapsNonUniqueCount;
	// multiple genes are mapped to the read, so it wouldn't be in the DGE.
	private int countMultiGeneMappingCount;

	// what genes receive unique mappings (not the original gene)
	private ObjectCounter<String> uniqueMapOtherGene;
	// what genes receive non-unique mappings (not the original gene)
	private ObjectCounter<String> nonUniqueMapOtherGene;
	// what contigs (other than the original) are mappings going to
	private ObjectCounter<String> otherContigs;
	private String noGeneTag;

	public GeneResult (final String originalGene, final String originalContig, final String noGeneTag) {
		this.originalGene=originalGene;
		this.originalContig=originalContig;
		this.countOriginalReads=0;
		this.countSameMapping=0;
		this.countNoGene=countNoGene;
		this.countDifferentUnique=0;
		this.countMapsNonUniqueCount=0;
		this.countMultiGeneMappingCount=0;
		this.uniqueMapOtherGene=new ObjectCounter<>();
		this.nonUniqueMapOtherGene=new ObjectCounter<>();
		this.otherContigs = new ObjectCounter<>();
		this.noGeneTag=noGeneTag;
	}

	/**
	 *
	 * @param genes The genes for each read.
	 */
	public void addMapping (final Collection <String> genes, final Collection<String> contigs, final int numReads) {
		if (this.originalGene.equals("AC002553.1"))
			System.out.println("STOP");

		boolean hasOriginalGene=genes.contains(this.originalGene);
		countOriginalReads++;

		// handle mappings to other contigs.
		contigs.remove(this.originalContig);
		for (String c: contigs)
			this.otherContigs.increment(c);

		// the only mapping is to the same gene.
		if (hasOriginalGene & numReads==1) {
			countSameMapping++;
			return;
		}
		// there's only one mapping to the original gene, but there are multiple reads.
		if (hasOriginalGene & genes.size()==1 & numReads>1) {
			countMapsNonUniqueCount++;
			return;
		}
		// maps to some new gene uniquely.
		if (!hasOriginalGene & genes.size()==1) {
			String gene = genes.iterator().next();
			if (gene.equals(this.noGeneTag)) {
				this.countNoGene++;
				return;
			} else {
				// otherwise you're a different gene.
				countDifferentUnique++;
				uniqueMapOtherGene.increment(genes.iterator().next());
				return;
			}
		}
		// there's more than 1 gene mapped.
		if (genes.size()>1) {
			countMultiGeneMappingCount++;
			genes.remove(this.originalGene);
			genes.remove(this.noGeneTag);
			for (String gene: genes)
				this.nonUniqueMapOtherGene.increment(gene);
			return;
		}
		throw new InvalidRequestStateException("Missed a case!");
	}

	public String getOriginalGene() {
		return originalGene;
	}

	public String getOriginalContig() {
		return originalContig;
	}

	public int getCountOriginalReads() {
		return countOriginalReads;
	}

	public int getCountSameMapping() {
		return countSameMapping;
	}

	public int getCountDifferentUnique() {
		return countDifferentUnique;
	}

	public int getCountMapsNonUniqueCount() {
		return countMapsNonUniqueCount;
	}

	public int getCountMultiGeneMappingCount() {
		return countMultiGeneMappingCount;
	}


	/**
	 * Get the unique mappings to genes that were not the original gene.
	 * @return
	 */
	public ObjectCounter<String> getUniqueMapOtherGene() {
		return uniqueMapOtherGene;
	}

	/**
	 * Get the non-unique mappings to genes that were not the original gene.
	 * @return
	 */
	public ObjectCounter<String> getNonUniqueMapOtherGene() {
		return nonUniqueMapOtherGene;
	}


	/**
	 * Get the non-unique mappings to genes that were not the original gene.
	 * @return
	 */
	public ObjectCounter<String> getOtherContigs() {
		return otherContigs;
	}


	@Override
	public String toString () {
		return "Gene [" + this.originalGene +"] contig [" + getOriginalContig() +"] original read count  [" + countOriginalReads+ "] same mapping [" + this.countSameMapping+"] different unique count [" + countDifferentUnique+"] non-unique read count [" + getCountMapsNonUniqueCount() +"] multimap gene counts [" + getCountMultiGeneMappingCount()
				+ "] other unique genes " + StringUtil.join(",", getUniqueMapOtherGene()) + " other non-unique genes " + StringUtil.join(",", getNonUniqueMapOtherGene()) + " other contigs " + StringUtil.join(",", getOtherContigs()) ;
	}

}
