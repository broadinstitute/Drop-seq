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

package org.broadinstitute.dropseqrna.eqtl;

import java.io.File;
import java.text.DecimalFormat;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.dropseqrna.annotation.GeneAnnotationReader;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.tools.DGEMatrix;
import org.broadinstitute.dropseqrna.eqtl.EqtlCovariate;
import org.broadinstitute.dropseqrna.eqtl.PrepareEqtlCovariates;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

import com.google.common.collect.Sets;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;
import picard.annotation.Gene;

public class CalculateXReactivationCovariate {

	private final Log log = Log.getInstance(CalculateXReactivationCovariate.class);
	private DecimalFormat decimalFormat = new DecimalFormat("#.0000");
	public static String FRACTION_X_COVAR_NAME="FRACTION_X";
	
	/**
	 * Calculate the X reactivation phenotype for each donor in the metaCellFile.
	 * This is the fraction of expression on the X chromosome compared to the fraction on the autosome. 
	 * 
	 * @param metaCellFile A file of meta cell expression.  Each row is a gene, each column a donor.
	 * @param annotationsFile A GTF or RefFlat file.
	 * @param contigGroupFile A YAML file containing the annotation groups each contig belongs to.  The file has a list of contig names, each of which has a list of annotation groups that contig belongs to.
	 * @param sequenceDictionaryFile A sequence dictionary file containing the expected contigs in the GTF.
	 * @param escapeGenesFile A list of genes that escape X inactivation.  This is optional, set to null to ignore.
	 * @return
	 */
	public EqtlCovariate getXReactivationCovariate(File metaCellFile, File annotationsFile, File contigGroupFile,
                                                   File sequenceDictionaryFile, File escapeGenesFile,
                                                   ValidationStringency validationStringency) {
		Set<String> getEscapeGenes = getEscapeGenes(escapeGenesFile);
		
		// load up the gene annotations.
		SamReader in = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE).open(sequenceDictionaryFile);
		SAMSequenceDictionary sd = in.getFileHeader().getSequenceDictionary();
		CloserUtil.close(in);
		OverlapDetector<Gene>  genes = GeneAnnotationReader.loadAnnotationsFile(annotationsFile, sd, validationStringency);
						
		// map the gene name to the contig.
		Map<String,String> geneContigMap = getGeneContigMap(genes);
		
		// get the contigs that belong to autosome or X groups.
		Map<ContigGroup, Set<String>> contigGroups=getContigGroups(contigGroupFile);
		Set<String> autosomeContigs = contigGroups.get(ContigGroup.AUTOSOME);
		Set<String> xContigs = contigGroups.get(ContigGroup.X);
		
		// load up the expression data.		
		DGEMatrix m = DGEMatrix.parseDenseFile(metaCellFile, "");
		
		// iterate over genes and gather up umi counts on autosomes vs X chromosome.
		ObjectCounter<String> autosomeCounts = new ObjectCounter<>();
		ObjectCounter<String> xCounts = new ObjectCounter<>();
		
		List<String> donorList=m.getCellBarcodes();
		
		for (String gene: m.getGenes()) {
			// ignore escape genes.
			if (getEscapeGenes.contains(gene)) continue;
			
			// is the gene on the autosomes or the X?
			String contig = geneContigMap.get(gene);
			if (contig==null) {
				log.warn("Gene " + gene + " in expression data not found in annotations - skipping!");
				continue;
			}
			
			double [] expression = m.getExpression(gene);
			// get the contig group and update.
			if (autosomeContigs.contains(contig))
				autosomeCounts=update(autosomeCounts, donorList, expression);
			if (xContigs.contains(contig))
				xCounts=update(xCounts, donorList, expression);							
		}
		
		// calculate and get final result.
		EqtlCovariate result = PrepareEqtlCovariates.generateEmptyCovariatesFile(metaCellFile);
		result = getFractionUmisX(result, autosomeCounts, xCounts);		
		return result;
	}
	
	private EqtlCovariate getFractionUmisX (EqtlCovariate covars, ObjectCounter<String> autosomeCounts, ObjectCounter<String> xCounts) {
		List<String> donorList=covars.donorNames();
		String [] fractionX = new String [donorList.size()];
		for (int i=0; i<donorList.size(); i++) {
			String donor=donorList.get(i);
			int x = xCounts.getCountForKey(donor);
			int a = autosomeCounts.getCountForKey(donor);
			double frac = (double)x / (double) a;
			fractionX[i] = decimalFormat.format(frac);			
		}
		covars.setValues(FRACTION_X_COVAR_NAME, fractionX);
		return covars;		
	}
	
	/**
	 * MetaData is in integer counts, so we Math.round the value from double->int.
	 * @param umiCounts
	 * @param donorList
	 * @param values
	 * @return
	 */
	private ObjectCounter<String> update (ObjectCounter<String> umiCounts, List<String> donorList, double [] values) {
		for (int i=0; i<donorList.size(); i++) {
			umiCounts.incrementByCount(donorList.get(i), (int) Math.round(values[i]));
		}
		return umiCounts;
	}
	
	private Map <ContigGroup, Set<String>> getContigGroups (File contigGroupFile) {
		// get the contig names for each group
		Map<String, Set<String>> contigGroups = ParseContigGroups.getContigGroupMapSimple(contigGroupFile);
		
		// validate expected groups are in the file.
		for (ContigGroup cg: ContigGroup.values()) {
			String cgName = cg.getName();
			if (!contigGroups.keySet().contains(cgName)) 
				log.error("The expected contig group ["+ ContigGroup.AUTOSOME.name() +"] not found in file " + contigGroupFile.getAbsolutePath());	
		}
		
		Set<String> autosomeContigs = contigGroups.get(ContigGroup.AUTOSOME.getName());
		Set<String> xContigs = contigGroups.get(ContigGroup.X.getName());
		
		// validate sets don't overlap
		Set<String> overlap = Sets.intersection(autosomeContigs, xContigs);
		if (overlap.size()>0) {			
			throw new IllegalArgumentException("Overlap in definition of autosomes and X chromosome contigs.  Contigs in both:"+overlap.toString());
		}
		Map <ContigGroup, Set<String>> result = new HashMap<>();
		result.put(ContigGroup.AUTOSOME, autosomeContigs);
		result.put(ContigGroup.X, xContigs);
		return (result);
	}

	private enum ContigGroup {
		AUTOSOME ("autosome"),
		X ("X");
		private final String name;
		
		ContigGroup (String name) {
			this.name=name;
		}
		
		public String getName() {
			return name;
		}
		
	}
	
	private Map<String,String> getGeneContigMap (OverlapDetector<Gene>  genes) {
		Map<String,String> result = new HashMap<>();
		
		for (Gene g: genes.getAll()) {
			result.put(g.getName(), g.getContig());
		}
		return result;
	}
	
	private Set<String> getEscapeGenes (File escapeGenes) {
		if (escapeGenes==null) return Collections.emptySet();
		return new HashSet<String> (ParseBarcodeFile.readCellBarcodeFile(escapeGenes));
	}
}
