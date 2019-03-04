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

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.PeekableIterator;

import java.util.*;

/**
 * Hold the DigitalAlleleCount results for many cells on a single SNP/Gene.
 * @author nemesh
 *
 */
public class MultiCellDigitalAlleleCounts {

	private final static String META_ANALYSIS_CELL="ALL_CELLS";

	private Map<String, DigitalAlleleCounts> dacMap;
	private final Interval snp;
	private final String gene;
	private boolean collapsed=false;

	/**
	 * Construct an object to hold the results for many cells on a single SNP/Gene.
	 * @param snp The interval of the SNP
	 * @param gene The name of the gene
	 */
	public MultiCellDigitalAlleleCounts(final Interval snp, final String gene) {
		this.snp = snp;
		this.gene=gene;
		dacMap = new HashMap<String, DigitalAlleleCounts>();
	}

	/**
	 * Collapses UMIs of all contained DigitalAlleleCount objects.
	 * After this operation is performed, no more data can be added to this object.
	 */
	public void collapseDACs (final int editDistance) {
		for (String key: dacMap.keySet()) {
			DigitalAlleleCounts dac = dacMap.get(key);
			dac = dac.collapseUMIs(editDistance);
			this.dacMap.put(key, dac);
		}
	}


	/**
	 * Only keep UMIs where the purity of the UMI is at least this threshold.
	 * @param threshold
	 */
	public void filterDataByUMIPurity (final double threshold) {
		for (String k: this.dacMap.keySet()) {
			DigitalAlleleCounts dac = this.dacMap.get(k);
			dac.filterDataByUMIPurity(threshold);
		}
	}

	public void add(final DigitalAlleleCounts dac) {
		if (collapsed)
			throw new IllegalArgumentException("Can't add more data to this object after UMI collapse!");
		if (dacMap.containsKey(dac.getCell()))
			throw new IllegalArgumentException("Already have cell [" + dac.getCell() +"] added to this collection!");
		dacMap.put(dac.getCell(), dac);
	}

	public Interval getSNP () {
		return this.snp;
	}

	public String getGene () {
		return this.gene;
	}

	/**
	 * Get a list of the cell names in this collection.
	 * These cells are sorted in alphabetical order.
	 * @return A list of cells
	 */
	public List<String> getCells () {
		List<String> cellIDs=new ArrayList<String>(dacMap.keySet());
		Collections.sort(cellIDs);
		return (cellIDs);
	}

	/**
	 * Get a DigitalAlleleCounts object for a particular cell
	 * @param cell The cell name to retrieve
	 * @return a DigitalAlleleCounts for a cell.
	 */
	public DigitalAlleleCounts getDigitalAlleleCounts(final String cell) {
		return dacMap.get(cell);
	}

	/**
	 * Constructs a meta cell containing all of the umi data across all individual cells.
	 * Allows you to calculate population wide statistics on all cells across a single SNP/gene.
	 * This changes UMI names to be a combination of the cell name and the UMI sequence.
	 * This is done because the same UMI sequence can exist in two different cells, and umi IDs are only unique within a cell.
	 * @return A new DigitalAlleleCounts object that contains the data from all cells added to this collection.
	 * If no data was added to the MultiCellDigitalAlleleCounts object, this will return null.
	 */
	public DigitalAlleleCounts getMetaAnalysis () {
		Iterator<DigitalAlleleCounts> iter =  dacMap.values().iterator();
		if (!iter.hasNext())
			return null;
		DigitalAlleleCounts first = iter.next();
		DigitalAlleleCounts result = new DigitalAlleleCounts(first.getSnpInterval(), first.getGene(), META_ANALYSIS_CELL, first.getBaseQualityThreshold());

		result=addDAC(result, first);
		while (iter.hasNext()) {
			DigitalAlleleCounts next=iter.next();
			result=addDAC(result, next);
		}
		return (result);
	}

	/**
	 * Constructs a meta cell containing all of the umi data across all individual cells in the provided list.
	 * Allows you to calculate population wide statistics on all cells across a single SNP/gene.
	 * This changes UMI names to be a combination of the cell name and the UMI sequence.
	 * This is done because the same UMI sequence can exist in two different cells, and umi IDs are only unique within a cell.
	 * @return A new DigitalAlleleCounts object that contains the data from all cells added to this collection.
	 * If no data was added to the MultiCellDigitalAlleleCounts object, this will return null.
	 */
	public DigitalAlleleCounts getMetaAnalysis (final String clusterID, final Collection<String> cells) {
		Iterator<DigitalAlleleCounts> iter =  dacMap.values().iterator();
		// short circuit.
		if (!iter.hasNext())
			return null;

		PeekableIterator<DigitalAlleleCounts> pIter= new PeekableIterator<DigitalAlleleCounts>(iter);
		DigitalAlleleCounts first = pIter.peek();
		DigitalAlleleCounts result = new DigitalAlleleCounts(first.getSnpInterval(), first.getGene(), clusterID, first.getBaseQualityThreshold());
		while (pIter.hasNext()) {
			DigitalAlleleCounts next=pIter.next();
			if (cells.contains(next.getCell()))
				result=addDAC(result, next);
		}
		pIter.close();
		return (result);
	}


	/**
	 * Adds data from a child DAC to a parent DAC and returns the parent.
	 * UMI names of the child are modified to have the child's cell ID.
	 * @param parent A parent DAC to add data to
	 * @param child A child DAC to get data from
	 * @return The parent DAC with additional data.  This modifies the parent variable.
	 */

	private DigitalAlleleCounts addDAC (final DigitalAlleleCounts parent, final DigitalAlleleCounts child) {
		for (String umi: child.umis()) {
			String key = getUMIKey(child.getCell(), umi);
			parent.addReadsForUMI(key, child.getReadsPerUMI(umi));
		}
		return (parent);
	}

	private String getUMIKey (final String cell, final String umi) {
		String key = cell + ":" + umi;
		return (key);
	}

}
