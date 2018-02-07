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


	// These reads map uniquely to a different gene and would show up in a DGE under a different gene.
	private int countDifferentUniqueGene;


	// LOSS OF EXPRESSION CATEGORIES

	// this read maps to one gene that isn't the original, and there are multiple reads so the new mapping would not be counted.
	private int countDifferentGeneNonUniqueCount;

	// these reads used to map uniquely to a gene, now they map to an intron/intergenic.
	private int countIntronicOrIntergenic;

	// the original gene mapping is present but read also maps to some other non-gene region, making the map quality of the gene read low so it wouldn't be in the DGE.
	// this read contains the original gene, but now there are multiple mappings.
	private int countSameGeneMapsNonUniqueCount;

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
		this.countIntronicOrIntergenic=0;
		this.countDifferentUniqueGene=0;
		this.countSameGeneMapsNonUniqueCount=0;
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
		boolean hasOriginalGene=genes.contains(this.originalGene);
		// check for the no-gene tag, set the flag, and remove it from the list.
		boolean hasReadNotOnAnyGene = genes.contains(this.noGeneTag);
		genes.remove(this.noGeneTag);
		// always incremented.
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

		// need to handle cases where a read maps to the same gene multiple times at low map quality.
		// kinda weird...I'm gonna call these intronic/intergenic, but they aren't really.
		if (hasOriginalGene & genes.size()==1 & numReads>1 & !hasReadNotOnAnyGene) {
			countSameGeneMapsNonUniqueCount++;
			return;
		}

		// there's only one mapping to the original gene, but there are multiple reads that map to no other genes.
		if (hasOriginalGene & genes.size()==1 & numReads>1 & hasReadNotOnAnyGene) {
			countSameGeneMapsNonUniqueCount++;
			return;
		}

		// doesn't map to a gene anymore.
		if (!hasOriginalGene & genes.size()==0 & hasReadNotOnAnyGene) {
			this.countIntronicOrIntergenic++;
			return;
		}
		// doesn't map to original gene but maps to some other gene.
		if (!hasOriginalGene & genes.size()==1)
			// maps to one other gene uniquely.
			if (numReads==1) {
				countDifferentUniqueGene++;
				uniqueMapOtherGene.increment(genes.iterator().next());
				return;
			} else { // maps to one other gene non-uniquely.
				countDifferentGeneNonUniqueCount++;
				return;
			}
		// there's more than 1 gene mapped.
		if (genes.size()>1) {
			genes.remove(this.originalGene);
			countMultiGeneMappingCount++;
			for (String gene: genes)
				this.nonUniqueMapOtherGene.increment(gene);
			return;
		}
		throw new IllegalStateException("Missed a case!");
	}

	/**
	 * The name of the gene
	 * @return
	 */
	public String getOriginalGene() {
		return originalGene;
	}

	/**
	 * The contig the gene was originally placed on.
	 * @return
	 */
	public String getOriginalContig() {
		return originalContig;
	}

	/**
	 * The number of original reads for this gene.
	 * @return
	 */
	public int getCountOriginalReads() {
		return countOriginalReads;
	}

	/**
	 * These reads map uniquely to the same read
	 * @return
	 */
	public int getCountSameMapping() {
		return countSameMapping;
	}

	/**
	 * These reads map uniquely to a different gene than the original read.
	 * @return
	 */
	public int getCountDifferentUniqueGene() {
		return countDifferentUniqueGene;
	}

	public int getCountDifferentGeneNonUniqueCount() {
		return countDifferentGeneNonUniqueCount;
	}

	/**
	 * The count of read that now map ONLY to a non-exonic part of the genome.
	 * @return
	 */
	public int getCountIntronicOrIntergenic() {
		return countIntronicOrIntergenic;
	}

	/**
	 * The count of reads that map to BOTH a gene and some other non-exonic region.
	 * @return
	 */
	public int getCountSameGeneMapsNonUniqueCount() {
		return countSameGeneMapsNonUniqueCount;
	}

	/**
	 * The count of reads that map to multiple gene exons.
	 * @return
	 */
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
		return "Gene [" + this.originalGene +"] contig [" + getOriginalContig() +"] original read count  [" + countOriginalReads+ "] same mapping [" + this.countSameMapping+"] different unique count [" + countDifferentUniqueGene+"] non-unique read count [" + getCountSameGeneMapsNonUniqueCount() +"] multimap gene counts [" + getCountMultiGeneMappingCount()
				+ "] other unique genes " + StringUtil.join(",", getUniqueMapOtherGene()) + " other non-unique genes " + StringUtil.join(",", getNonUniqueMapOtherGene()) + " other contigs " + StringUtil.join(",", getOtherContigs()) ;
	}

}
