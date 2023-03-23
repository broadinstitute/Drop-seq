package org.broadinstitute.dropseqrna.barnyard;

import java.util.Collections;
import java.util.Set;

import org.broadinstitute.dropseqrna.utils.ObjectCounter;

public class ChimericUmi {
	
	private final String molecularBarcode;	
	private final Set<String> chimericGenes;
	
	public enum CHIMERIC_STRATEGY {
    	REMOVE_ALL, RETAIN_MOST_SUPPORTED;
    }

	/**
	 * Construct a ChimericUmi object for a single cell.
	 * This object receives the number of reads supporting each gene for a cell and UMI sequence, and based on the strategy
	 * remembers the set of genes that are chimeric.
	 * @param molecularBarcode The molecular barcode
	 * @param genes A collection of gene symbols with the number of supporting reads for this cell/molecular barcode.
	 * @param strategy The strategy to discover chimeric reads.  REMOVE_ALL removes all genes for this UMI if there is more than 1 gene.  RETAIN_MOST_SUPPORTED retains
	 * the gene with the highest read count and flags all other genes as chimeric.  If the most supported gene is ambiguous (there are multiple genes with the highest count) then 
	 * all genes are chimeric. 
	 */
	public ChimericUmi(String molecularBarcode, ObjectCounter<String> genes, CHIMERIC_STRATEGY strategy) {
		this.molecularBarcode=molecularBarcode;
		this.chimericGenes=getChimericGenes(genes,strategy);
	}

	public String getMolecularBarcode() {
		return molecularBarcode;
	}
	
	public Set<String> getChimericGenes() {
		return chimericGenes;
	}
	
	/**
	 * Find the set of 0 or more genes where the UMI 
	 * @return
	 */
	static Set<String> getChimericGenes (ObjectCounter<String> genes, CHIMERIC_STRATEGY strategy) {    		
		switch (strategy) {
		case REMOVE_ALL: return (getChimericGenesRemoveAll(genes));
		case RETAIN_MOST_SUPPORTED: return(getChimericGenesRetainMostSupported(genes));
		default: 
			throw new IllegalArgumentException("Chimeric Strategy not supported");
		}			
	}
	
	/**
	 * Simple strategy to filter chimeric genes.  If multiple genes share the same UMI, all genes are chimeric.
	 * @param genes
	 * @return
	 */
	static Set <String> getChimericGenesRemoveAll(ObjectCounter<String> genes) {
		if (genes.getSize()<2) return Collections.emptySet();
		return (Set<String>) genes.getKeys();
	}
	
	/**
	 * If there are multiple genes with varying degrees of support, flag all but the most
	 * supported gene as chimeric.  If multiple genes are the most supported, the result is ambiguous
	 * and all genes for this UMI are considered chimeric.
	 * @param genes A set of genes with read counts
	 * @return A set of chimeric gene symbols.  This set can be empty.
	 */
	static Set<String> getChimericGenesRetainMostSupported(ObjectCounter<String> genes) {
		// no chimeric genes.
		if (genes.getSize()<2) return Collections.emptySet();
		// at least one chimeric gene.
		String mostSupportedGene=genes.getMax();
		int counts = genes.getCountForKey(mostSupportedGene);
		int genesAtSupport=genes.getNumberOfSize(counts);
		if (genesAtSupport==1) {
			Set<String> chimeric = (Set<String>) genes.getKeys();
			chimeric.remove(mostSupportedGene);
			return chimeric;
		}			
		// ambiguous, return all genes as chimeric.
		return (Set<String>)genes.getKeys();		
	}

	
}