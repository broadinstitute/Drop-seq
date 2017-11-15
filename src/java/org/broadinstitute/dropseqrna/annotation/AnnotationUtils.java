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
package org.broadinstitute.dropseqrna.annotation;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import picard.annotation.Gene;
import picard.annotation.Gene.Transcript;
import picard.annotation.Gene.Transcript.Exon;
import picard.annotation.LocusFunction;

public class AnnotationUtils {

	@SuppressWarnings("unused")
	private final Log log = Log.getInstance(AnnotationUtils.class);
	
	private static AnnotationUtils singleton=null;
	
	private static Map<LocusFunction, Integer> functionScores;
	
	private AnnotationUtils() {
		functionScores = new HashMap<LocusFunction, Integer>();
		functionScores.put(LocusFunction.CODING, new Integer (4));
		functionScores.put(LocusFunction.UTR, new Integer (3));
		functionScores.put(LocusFunction.INTRONIC, new Integer (2));
		functionScores.put(LocusFunction.INTERGENIC, new Integer (1));
	}
	
	public static AnnotationUtils getInstance() {
		if (singleton==null) {
			singleton=new AnnotationUtils();
		}
		return singleton;
	}
	
	/**
	 * Given a BAM record and Gene object, output the functional class of the read.
	 * When a read overlaps multiple locus functions, the top function is selected in the order 1) exon 2) intron 3) UTR 4) intergenic.
	 * This means a read that overlaps both an exon and an intron is annotated as exonic
	 * @param rec The sam Record
	 * @param geneOverlapDetector An overlap detector, such as that produced by GeneAnnotationReader.
	 * @return The function of the read.
	 */
	public LocusFunction getLocusFunctionForRead (SAMRecord rec, OverlapDetector<Gene> geneOverlapDetector) {
		final Interval readInterval = new Interval(rec.getReferenceName(), rec.getAlignmentStart(), rec.getAlignmentEnd());
		final Collection<Gene> overlappingGenes = geneOverlapDetector.getOverlaps(readInterval);
		
		
        final List<AlignmentBlock> alignmentBlocks = rec.getAlignmentBlocks();
        
        LocusFunction [] blockSummaryFunction = new LocusFunction[alignmentBlocks.size()];
        
        for (int i=0; i<alignmentBlocks.size(); i++) {
        	AlignmentBlock alignmentBlock =alignmentBlocks.get(i);
        	LocusFunction [] blockFunctions=getLocusFunctionsByBlock(alignmentBlock, overlappingGenes);
        	LocusFunction blockFunction = getLocusFunction(blockFunctions);
        	blockSummaryFunction[i]=blockFunction;
        }
        
        LocusFunction readFunction = getLocusFunction(blockSummaryFunction);
		return readFunction;
	}
	
	
	/**
	 * Get the LocusFunction on a gene by gene basis for the alignment blocks of this read.  
	 * If a read is split into multiple blocks, each block should point to the same gene - if blocks point to 
	 * different genes (which can't happen, you can't splice different genes together), then that gene result should not be returned.
	 * Instead, return null in those cases. 
	 * @param alignmentBlocks
	 * @param g
	 * @return
	 */
	private LocusFunction getLocusFunctionForRead (SAMRecord rec, Gene g) {
		List<AlignmentBlock> alignmentBlocks = rec.getAlignmentBlocks();
		
        LocusFunction [] blockSummaryFunction = new LocusFunction[alignmentBlocks.size()];
        Set<Gene> temp = new HashSet<Gene>();
        temp.add(g);
        
        for (int i=0; i<alignmentBlocks.size(); i++) {
        	AlignmentBlock alignmentBlock =alignmentBlocks.get(i);
        	
        	LocusFunction [] blockFunctions=getLocusFunctionsByBlock(alignmentBlock, temp);
        	LocusFunction blockFunction = getLocusFunction(blockFunctions);
        	blockSummaryFunction[i]=blockFunction;
        }
        
        LocusFunction readFunction = getLocusFunction(blockSummaryFunction);
		return readFunction;
	}
	
	
		
	public Map<Gene, LocusFunction> getLocusFunctionForReadByGene (SAMRecord rec, OverlapDetector<Gene> geneOverlapDetector) {
		Map<Gene, LocusFunction> result = new HashMap<Gene, LocusFunction>();
		final Interval readInterval = new Interval(rec.getReferenceName(), rec.getAlignmentStart(), rec.getAlignmentEnd());
		final Collection<Gene> overlappingGenes = geneOverlapDetector.getOverlaps(readInterval);
		
        for (Gene g: overlappingGenes) {
        	LocusFunction f = getLocusFunctionForRead(rec, g);
            result.put(g, f);
        }
        	
        
        return result;
	}
	
	
	/**
	 * For each alignment block, see which genes have exons that intersect it.
	 * Only count genes where the alignment blocks are consistent - either the alignment blocks point to the same gene, or 
	 * an alignment block points to a gene and the other to no genes/exons in the set. 
	 * 
	 * @param rec
	 * @param genes A set of genes that the read originally overlaps.
	 * @param allowMultipleGenes.  If false, and a read overlaps multiple gene exons, then none of the genes are returned.  If true, return the set of all genes.
	 * @return
	 */
	public Set<Gene> getConsistentExons (SAMRecord rec, Set<Gene> genes, boolean allowMultiGeneReads) {
		
		Set<Gene> result = new HashSet<Gene>();
		String refName = rec.getReferenceName();
		List<AlignmentBlock> alignmentBlocks = rec.getAlignmentBlocks();
		for (AlignmentBlock b: alignmentBlocks) {
			Set<Gene> blockGenes = getAlignmentBlockonGeneExon(refName, b, genes);
			// if result is not empty and blockGenes isn't empty, intersect the set and set the result as this new set.
			if (result.size()>0 && blockGenes.size()>0) {
				if (allowMultiGeneReads) result.addAll(blockGenes);
				if (!allowMultiGeneReads)  result.retainAll(blockGenes);				
			} else {
				// if blockGenes is populated and you're here, then result is empty, so set result to these results
				result=blockGenes;
			}
		}
		return (result);		
	}
	
	private Set<Gene> getAlignmentBlockonGeneExon(String refName, AlignmentBlock b, Set<Gene> genes) {
		Set<Gene> result = new HashSet<Gene>();
		for (Gene g: genes) {
			if (getAlignmentBlockOverlapsExon(refName, b, g)) {
				result.add(g);
			}
		}
		return (result);
	}
	
	private boolean getAlignmentBlockOverlapsExon (String refName, AlignmentBlock b, Gene g) {
		
		Interval ib = getInterval(refName, b);
		
		for (Transcript t: g) {
			for (Exon e: t.exons) {
				Interval ei = getInterval(refName, e);
				if (ib.intersects(ei)) {
					return true;
				}
			}
		}
		return false;
	}
	
	private Interval getInterval (String refName, Exon e) {
		Interval i = new Interval(refName, e.start,e.end);
		return (i);
	}
	
	private Interval getInterval (String refName, AlignmentBlock b) {
		int s = b.getReferenceStart();
		int e = s + b.getLength()-1;
		Interval i = new Interval(refName, s,e);
		return (i);
	}
	
	
	
	public LocusFunction getLocusFunction (Collection<LocusFunction> locusFunctions) {
		if (locusFunctions.size()==0) return LocusFunction.INTERGENIC;
		LocusFunction[] array = locusFunctions.toArray(new LocusFunction[locusFunctions.size()]);
		return getLocusFunction(array);
	}
	
	public LocusFunction getLocusFunction (LocusFunction[] locusFunctions) {
		if (locusFunctions==null || locusFunctions.length==0) return (null);
		
		int bestScore=Integer.MIN_VALUE;
		LocusFunction bestFunction=LocusFunction.INTERGENIC;
		for (LocusFunction f: locusFunctions) {
			int score = functionScores.get(f);
			if (score>bestScore) {
				bestScore=score;
				bestFunction=f;
			}
		}
		return bestFunction;
	}
	
	
	public LocusFunction [] getLocusFunctionsByBlock (AlignmentBlock b, Collection<Gene> overlappingGenes) {
		// Get functional class for each position in the alignment block.
        final LocusFunction[] locusFunctions = new LocusFunction[b.getLength()];

        // By default, if base does not overlap with rRNA or gene, it is intergenic.
        Arrays.fill(locusFunctions, 0, locusFunctions.length, LocusFunction.INTERGENIC);

        for (final Gene gene : overlappingGenes) {
            for (final Gene.Transcript transcript : gene) {
                transcript.assignLocusFunctionForRange(b.getReferenceStart(), locusFunctions);
            }
        }
        return locusFunctions;
	}
	
	public LocusFunction[] getLocusFunctionsByInterval (Interval i, Collection<Gene> genes) {
		// Get functional class for each position in the alignment block.
        final LocusFunction[] locusFunctions = new LocusFunction[i.length()];

        // By default, if base does not overlap with rRNA or gene, it is intergenic.
        Arrays.fill(locusFunctions, 0, locusFunctions.length, LocusFunction.INTERGENIC);

        for (final Gene gene : genes) {
            for (final Gene.Transcript transcript : gene) {
                transcript.assignLocusFunctionForRange(i.getStart(), locusFunctions);
            }
        }
        return locusFunctions;
	}
	
	
	Map<String, String> parseOptionalFields(String optional) {
		Map<String, String> result = new HashMap<String, String>();
		String [] o = optional.split(";");
		for (String s: o) {
			s=s.replaceAll("\"", "");
			s=s.trim();
			String [] z= s.split(" ");
			String k = z[0];
			String v = z[1];
			result.put(k, v);
		}
		return (result);
	}
	
	public static String strandToString(boolean isPositiveStrand) {
		if (isPositiveStrand) return "+";
		return "-";
	}
}
