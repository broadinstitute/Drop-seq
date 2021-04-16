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

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.broadinstitute.dropseqrna.utils.ObjectCounter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import picard.annotation.LocusFunction;

/**
 * A pileup of bases and qualities at a SNP in the context of a UMI (gene/cell/molecular barcode).
 * @author nemesh
 *
 */
public class SNPUMIBasePileup extends SNPBasePileUp {

	private static final Log log = Log.getInstance(SNPUMIBasePileup.class);

	private final String gene;
	private final String cell;
	private final String molecularBarcode;
	private Set<LocusFunction> locusFunctionSet;

	public SNPUMIBasePileup (final Interval snpInterval, final String gene, final String cell, final String molecularBarcode) {
		super(snpInterval);
		this.gene=gene;
		this.cell=cell;
		this.molecularBarcode=molecularBarcode;
		this.locusFunctionSet=new HashSet<>();
	}

	public String getGene() {
		return gene;
	}

	public String getCell() {
		return cell;
	}

	public String getMolecularBarcode() {
		return molecularBarcode;
	}


	public void addLocusFunction (final LocusFunction f) {
		this.locusFunctionSet.add(f);
	}


	/**
	 * For each read, extract the base and base quality at the snp position of this object and add to the pileup.
	 * @param r The SAMRecord to extract the base and quality for this SNP position.
	 */
	@Override
	public void addRead (final SAMRecord r) {
		byte [] result = getBaseAndQualityOverlappingInterval(r);
		if (result.length!=0)
			addBaseAndQuality(result[0], result[1]);
	}

	public Set<LocusFunction> getLocusFunctions() {
		return this.locusFunctionSet;
	}
	
	public double getUMIPurity () {
		ObjectCounter<Character> counts = new ObjectCounter<>();
		this.getBasesAsCharacters().stream().forEach(x -> counts.increment(x));		
		List<Character> keys = counts.getKeysOrderedByCount(true);
		int first = counts.getCountForKey(keys.get(0));
		int total = counts.getTotalCount();				
		double result = (double) first / (double) (total);
		return result;
	}

	@Override
	public String toString () {
		StringBuilder b = new StringBuilder();
		b.append("snp [" + this.getSNPInterval().toString() +"] ");
		b.append("gene [" + this.gene +"] ");
		b.append("cell [" + this.cell +"] ");
		b.append("molBC [" + this.molecularBarcode +"] ");
		b.append(this.getBasesAsCharacters().toString()+ " ");
		b.append(this.getQualities().toString()+" ");
		if (this.locusFunctionSet!=null && this.locusFunctionSet.size()>0)
			b.append("locus function " + this.locusFunctionSet.toString()+"");
		return b.toString();
	}

}
