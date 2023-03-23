package org.broadinstitute.dropseqrna.barnyard;

import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.dropseqrna.barnyard.ChimericUmi.CHIMERIC_STRATEGY;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

import htsjdk.samtools.util.Log;

/**
 * A collection of all chimeric UMIs / genes for a single cell.
 * @author nemesh
 *
 */
public class ChimericUmiCollection {

	private static final Log log = Log.getInstance(ChimericUmiCollection.class);
	
	private final String cellBarcode;
	// map from a UMI to the genes that are chimeric for that UMI.
	private Map<String, ChimericUmi> map;
	// key is the molecular barcode, value is the number of reads per gene.
	
	private Set<String> problematicUmis;

	private Map<String, ObjectCounter<String>> umiGeneSupport;
	private boolean isFinal = false;
	private final ChimericUmi.CHIMERIC_STRATEGY strategy;
	private int totalUMIs=0;
	
	
	public ChimericUmiCollection (String cellBarcode, ChimericUmi.CHIMERIC_STRATEGY strategy) {
		this.cellBarcode=cellBarcode;
		this.strategy=strategy;
		totalUMIs=0;
		map = new HashMap<>();
		// A map from the UMI sequence to the genes that belong to that UMI with the number of reads per gene.
		// This will be removed when the data is finalized to save memory.
		umiGeneSupport = new HashMap<String, ObjectCounter<String>>();
		this.problematicUmis=new HashSet<>();
	}
		
	/**
	 * Add data from a UMICollection directly.
	 * @param u
	 */
	public void add (UMICollection u) {
		add(u.getGeneName(), u.getMolecularBarcodeCounts());
	}
	
	/**
	 * Add data for a single gene and the collection of molecular barcodes for that gene with their corresponding read counts
	 * @param gene
	 * @param molBCCounts
	 */
	public void add (String gene, ObjectCounter<String> molBCCounts) {
		if (isFinal) 
			throw new IllegalStateException("This object has been finalized, no more data can be added");
		
		for (String molBC: molBCCounts.getKeys()) {
			
    		int support = molBCCounts.getCountForKey(molBC);
    		ObjectCounter<String> geneSupportCounts = umiGeneSupport.get(molBC);
    		if (geneSupportCounts==null) {
    			geneSupportCounts=new ObjectCounter<>();
    			umiGeneSupport.put(molBC, geneSupportCounts);
    		}
    		geneSupportCounts.incrementByCount(gene, support);
    		totalUMIs = getTotalUMIs() + 1;
    	}
	}
		
	/**
	 * Get a list of 0 or more genes that are chimeric for the requested UMI.
	 * When this method is called the object is finalized and no more data can be added. 
	 * @param umi A molecular barcode to query
	 * @return A set with 0 or more entries of chimeric genes.  If the UMI does not exist for this cell, return an empty set of genes.
	 */
	public Set<String> getChimericGenes (String umi) {		
		if (isFinal==false) 
			buildChimericUmis();
		ChimericUmi cu = map.get(umi);
		if (cu==null)
			return Collections.emptySet();
		return cu.getChimericGenes();
	}
	
	/**
	 * Checks if a UMI is chimeric, either via the chimerism test or via additionally registered UMIs that are problematic.
	 * @param umi The molecular barcode to test.
	 * @return true if the read is chimeric, or has been registered as problematic.
	 */
	public boolean isChimericOrProblematic (String umi, String gene) {
		if (isFinal==false) 
			buildChimericUmis();
		if (this.problematicUmis.contains(umi))
			return true;
		Set<String> chimericGenes = getChimericGenes(umi);
		return chimericGenes.contains(gene);		
	}
	
	/**
	 * Record that a UMI is problematic for this cell by a different test than the chimeric read test.
	 * @param umi
	 */
	public void registerProblematicUmi (String umi) {
		this.problematicUmis.add(umi);
	}
	
	public Set<String> getChimericUmis () {
		return map.keySet();
	}
	
	public Set<String> getProblematicUmis() {
		return problematicUmis;
	}

		
	
	/**
	 * Finalize the data set.  Chimeric UMIs are discovered from the input data.
	 * Input data is discarded to save memory, and no more data can be added.
	 * It is recommended to finalize each cell barcode's data as it is entered to save memory 
	 * since only the genes that are chimeric for a given molecular barcode are retained.
	 */	
	public void buildChimericUmis () {		
		for (String molBC: umiGeneSupport.keySet()) {
			ChimericUmi cu = new ChimericUmi(molBC, umiGeneSupport.get(molBC), this.strategy);
			map.put(molBC, cu);
		}
		this.isFinal=true;
		this.umiGeneSupport=null;
	}

	public String getCellBarcode() {
		return cellBarcode;
	}

	public int getTotalUMIs() {		
		return totalUMIs;
	}
	
	/**
	 * Find the number of chimeric UMI/gene pairs for this cell.
	 * @return
	 */
	public int getTotalUMisChimeric() {
		if (isFinal==false) 
			buildChimericUmis();
		int count=0;
		for (String umi: this.map.keySet()) {
			count+=map.get(umi).getChimericGenes().size();
		}
		return (count);
	}
	
	
	
}
