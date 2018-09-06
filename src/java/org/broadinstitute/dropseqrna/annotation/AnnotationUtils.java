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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.collections4.map.HashedMap;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;
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
		functionScores = new HashMap<>();
		functionScores.put(LocusFunction.CODING, new Integer (5));
		functionScores.put(LocusFunction.UTR, new Integer (4));
		functionScores.put(LocusFunction.INTRONIC, new Integer (3));
		functionScores.put(LocusFunction.RIBOSOMAL, new Integer (2));
		functionScores.put(LocusFunction.INTERGENIC, new Integer (1));
	}

	public static AnnotationUtils getInstance() {
		if (singleton==null)
			singleton=new AnnotationUtils();
		return singleton;
	}


	/**
	 * For a read, split into the alignment blocks.
	 * For each alignment block, determine the genes the block overlaps, and the locus function of those genes.
	 *
	 * Each alignment block can generate a different functional annotation on the same gene.  These should be retained.
	 * For example, block one can align to an exon, and block two to an intron, and both those annotations are retained.
	 *
	 * Only retain genes where alignment blocks all reference that gene.  If block one refers to genes A,B and block two to gene A only, then only retain gene A.
	 *
	 */
	public Map<Gene, List<LocusFunction>> getFunctionalDataForRead (final SAMRecord rec, final OverlapDetector<Gene> geneOverlapDetector) {
		List<AlignmentBlock> alignmentBlocks = rec.getAlignmentBlocks();

		Map<AlignmentBlock, Map<Gene, List<LocusFunction>>> map = new HashMap<>();
		// gather the locus functions for each alignment block.
		for (AlignmentBlock block: alignmentBlocks) {
			Interval interval = getInterval(rec.getReferenceName(), block);
			Map<Gene, List<LocusFunction>> locusFunctionsForGeneMap = getFunctionalDataForInterval(interval, geneOverlapDetector);
			map.put(block, locusFunctionsForGeneMap);
		}
		// simplify genes by only using genes that are common to all alignment blocks.
		Map<Gene, List<LocusFunction>> result = simplifyFunctionalDataAcrossAlignmentBlocks(map);
		return result;
	}

	/**
	 * Get a list of functional data - genes and their locus functions for an interval.
	 * A gene can appear multiple times in the output with different locus functions.
	 * @param interval The interval to interrogate
	 * @param geneOverlapDetector The gene model to interrogate.
	 * @return A map from each gene that is overlapped to its functional annotation(s).
	 */
	private Map<Gene, List<LocusFunction>> getFunctionalDataForInterval (final Interval interval, final OverlapDetector<Gene> geneOverlapDetector) {
		Map<Gene, List<LocusFunction>> result = new HashMap<>();

		final Collection<Gene> overlappingGenes = geneOverlapDetector.getOverlaps(interval);
		for (Gene g: overlappingGenes) {
			List<LocusFunction> locusFunctionsForGene = new ArrayList<>();

			LocusFunction [] locusFunctionArray = getLocusFunctionsByInterval(interval, g);
			// simplify to the unique list of functions.
			ObjectCounter<LocusFunction> o = new ObjectCounter<>();
			for (LocusFunction f: locusFunctionArray)
				o.increment(f);
			for (LocusFunction f: o.getKeys())
				locusFunctionsForGene.add(f);
			result.put(g,locusFunctionsForGene);
		}
		return (result);
	}

	/**
	 * Given the genes and functional annotations for each alignment block, simplify results to the read.
	 * Only retain genes that are common across alignment blocks.
	 * @param map
	 * @return
	 */
	private Map<Gene, List<LocusFunction>> simplifyFunctionalDataAcrossAlignmentBlocks (final Map<AlignmentBlock, Map<Gene, List<LocusFunction>>> map) {
		// if map is empty.
		if (map.size()==0)
			return Collections.emptyMap();
		// short circuit if there's only 1 alignment block.
		if (map.size()==1)
			return map.values().iterator().next();

		// initialize with the genes from the first alignment block.
		Iterator<Map<Gene, List<LocusFunction>>> iter = map.values().iterator();
		Set<Gene> commonGenes = iter.next().keySet();
		while (iter.hasNext()) {
			Set<Gene> next = iter.next().keySet();
			commonGenes.retainAll(next);
		}

		// walk through alignment blocks and retain genes in the common set.
		Map<Gene, List<LocusFunction>> result = new HashedMap<>();
		for (AlignmentBlock b: map.keySet())
			for (Gene g: commonGenes) {
				List<LocusFunction> tempResult = result.get(g);
				if (tempResult==null) tempResult = new ArrayList<>();
				List<LocusFunction> lf = map.get(b).get(g);
				tempResult.addAll(lf);
				result.put(g, tempResult);
			}

		// remove any repeats of a locus function for a gene - so if two blocks are both coding, only have 1 coding in the result.
		for (Gene g: result.keySet()) {
			List<LocusFunction> lf = result.get(g);
			Set<LocusFunction> s = new HashSet<>(lf);
			lf.clear();
			lf.addAll(s);
		}

		return result;
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
	private LocusFunction getLocusFunctionForRead (final SAMRecord rec, final Gene g) {
		List<AlignmentBlock> alignmentBlocks = rec.getAlignmentBlocks();

        LocusFunction [] blockSummaryFunction = new LocusFunction[alignmentBlocks.size()];
        Set<Gene> temp = new HashSet<>();
        temp.add(g);

        for (int i=0; i<alignmentBlocks.size(); i++) {
        	AlignmentBlock alignmentBlock =alignmentBlocks.get(i);

        	LocusFunction [] blockFunctions=getLocusFunctionsByBlock(alignmentBlock, temp);
        	LocusFunction blockFunction = getLocusFunction(blockFunctions, false);
        	blockSummaryFunction[i]=blockFunction;
        }

        LocusFunction readFunction = getLocusFunction(blockSummaryFunction, false);
		return readFunction;
	}



	public Map<Gene, LocusFunction> getLocusFunctionForReadByGene (final SAMRecord rec, final OverlapDetector<Gene> geneOverlapDetector) {
		Map<Gene, LocusFunction> result = new HashMap<>();
		final Interval readInterval = new Interval(rec.getReferenceName(), rec.getAlignmentStart(), rec.getAlignmentEnd(), rec.getReadNegativeStrandFlag(), rec.getReadName());
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
	 * Note that this is not strand specific, which may cause problems if a read overlaps two genes on opposite strands that have overlapping exons.
	 * @param rec
	 * @param genes A set of genes that the read originally overlaps.
	 * @param allowMultipleGenes.  If false, and a read overlaps multiple gene exons, then none of the genes are returned.  If true, return the set of all genes.
	 * @return
	 */
	//TODO: if this is used in the future, make it strand specific.
	public Set<Gene> getConsistentExons (final SAMRecord rec, final Set<Gene> genes, final boolean allowMultiGeneReads) {

		Set<Gene> result = new HashSet<>();
		String refName = rec.getReferenceName();
		List<AlignmentBlock> alignmentBlocks = rec.getAlignmentBlocks();
		for (AlignmentBlock b: alignmentBlocks) {
			Set<Gene> blockGenes = getAlignmentBlockonGeneExon(refName, b, genes);
			result.addAll(blockGenes);
			/*
			// if result is not empty and blockGenes isn't empty, intersect the set and set the result as this new set.
			if (result.size()>0 && blockGenes.size()>0) {
				if (allowMultiGeneReads) result.addAll(blockGenes);
				if (!allowMultiGeneReads)  result.retainAll(blockGenes);
			} else
				// if blockGenes is populated and you're here, then result is empty, so set result to these results
				result=blockGenes;
			*/
		}
		if (!allowMultiGeneReads & result.size()>1) return new HashSet<>();
		return result;
	}



	private Set<Gene> getAlignmentBlockonGeneExon(final String refName, final AlignmentBlock b, final Set<Gene> genes) {
		Set<Gene> result = new HashSet<>();
		for (Gene g: genes)
			if (getAlignmentBlockOverlapsExon(refName, b, g))
				result.add(g);
		return (result);
	}

	private boolean getAlignmentBlockOverlapsExon (final String refName, final AlignmentBlock b, final Gene g) {

		Interval ib = getInterval(refName, b);

		for (Transcript t: g)
			for (Exon e: t.exons) {
				Interval ei = getInterval(refName, e);
				if (ib.intersects(ei))
					return true;
			}
		return false;
	}

	private Interval getInterval (final String refName, final Exon e) {
		Interval i = new Interval(refName, e.start,e.end);
		return (i);
	}

	private Interval getInterval (final String refName, final AlignmentBlock b) {
		int s = b.getReferenceStart();
		int e = s + b.getLength()-1;
		Interval i = new Interval(refName, s,e);
		return (i);
	}



	public LocusFunction getLocusFunction (final Collection<LocusFunction> locusFunctions, final boolean conservative) {
		if (locusFunctions.size()==0) return LocusFunction.INTERGENIC;
		LocusFunction[] array = locusFunctions.toArray(new LocusFunction[locusFunctions.size()]);
		return getLocusFunction(array, conservative);
	}

	/**
	 * Summarize the locus functions that are at base-by-base level to a single annotation.
	 *
	 * @param locusFunctions
	 * @param conservative If true, only return a LocusFunction if all locusFunctions are the same, otherwise return null.  If false, return the "best" annotation, where annotations like coding are preferred over intronic,intergenic.
	 * @return
	 */
	public LocusFunction getLocusFunction (final LocusFunction[] locusFunctions, final boolean conservative) {
		if (locusFunctions==null || locusFunctions.length==0) return (null);
		if (conservative) {
			LocusFunction f = locusFunctions[0];
			for (LocusFunction lf: locusFunctions)
				if (lf!=f) return null;
			return f;
		}

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


	public LocusFunction [] getLocusFunctionsByBlock (final AlignmentBlock b, final Collection<Gene> overlappingGenes) {
		// Get functional class for each position in the alignment block.
        final LocusFunction[] locusFunctions = new LocusFunction[b.getLength()];

        // By default, if base does not overlap with rRNA or gene, it is intergenic.
        Arrays.fill(locusFunctions, 0, locusFunctions.length, LocusFunction.INTERGENIC);

        for (final Gene gene : overlappingGenes)
			for (final Gene.Transcript transcript : gene)
				transcript.assignLocusFunctionForRange(b.getReferenceStart(), locusFunctions);
        return locusFunctions;
	}

	public LocusFunction [] getLocusFunctionsByBlock (final Interval interval, final Collection<Gene> overlappingGenes) {
		// Get functional class for each position in the alignment block.
        final LocusFunction[] locusFunctions = new LocusFunction[interval.length()];

        // By default, if base does not overlap with rRNA or gene, it is intergenic.
        Arrays.fill(locusFunctions, 0, locusFunctions.length, LocusFunction.INTERGENIC);

        for (final Gene gene : overlappingGenes)
			for (final Gene.Transcript transcript : gene)
				transcript.assignLocusFunctionForRange(interval.getStart(), locusFunctions);
        return locusFunctions;
	}



	/**
	 * Generates the locus function at each base of the interval.
	 * @param interval
	 * @param overlappingGenes
	 * @return
	 */
	public LocusFunction[] getLocusFunctionsByInterval (final Interval i, final Gene gene) {
		// Get functional class for each position in the alignment block.
        final LocusFunction[] locusFunctions = new LocusFunction[i.length()];

        // By default, if base does not overlap with rRNA or gene, it is intergenic.
        Arrays.fill(locusFunctions, 0, locusFunctions.length, LocusFunction.INTERGENIC);

		for (final Gene.Transcript transcript : gene)
			transcript.assignLocusFunctionForRange(i.getStart(), locusFunctions);
        return locusFunctions;
	}


	Map<String, String> parseOptionalFields(final String optional) {
		Map<String, String> result = new HashMap<>();
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

	public static String strandToString(final boolean isPositiveStrand) {
		if (isPositiveStrand) return "+";
		return "-";
	}


	/**
	 * Given a BAM record and Gene object, output the functional class of the read.
	 * When a read overlaps multiple locus functions, the top function is selected in the order 1) exon 2) intron 3) UTR 4) intergenic.
	 * This means a read that overlaps both an exon and an intron is annotated as exonic
	 * @param rec The sam Record
	 * @param geneOverlapDetector An overlap detector, such as that produced by GeneAnnotationReader.
	 * @return The function of the read.
	 */
	/*
	public LocusFunction getLocusFunctionForRead (final SAMRecord rec, final OverlapDetector<Gene> geneOverlapDetector) {
		final Interval readInterval = new Interval(rec.getReferenceName(), rec.getAlignmentStart(), rec.getAlignmentEnd());
		final Collection<Gene> overlappingGenes = geneOverlapDetector.getOverlaps(readInterval);

        final List<AlignmentBlock> alignmentBlocks = rec.getAlignmentBlocks();

        LocusFunction [] blockSummaryFunction = new LocusFunction[alignmentBlocks.size()];

        for (int i=0; i<alignmentBlocks.size(); i++) {
        	AlignmentBlock alignmentBlock =alignmentBlocks.get(i);
        	LocusFunction [] blockFunctions=getLocusFunctionsByBlock(alignmentBlock, overlappingGenes);
        	LocusFunction blockFunction = getLocusFunction(blockFunctions, false);
        	blockSummaryFunction[i]=blockFunction;
        }

        LocusFunction readFunction = getLocusFunction(blockSummaryFunction, false);
		return readFunction;
	}
	*/



}
