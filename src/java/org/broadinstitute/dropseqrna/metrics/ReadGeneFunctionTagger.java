/*
 * MIT License
 *
 * Copyright 2026 Broad Institute
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
package org.broadinstitute.dropseqrna.metrics;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.dropseqrna.annotation.AnnotationUtils;
import org.broadinstitute.dropseqrna.barnyard.Utils;
import picard.annotation.Gene;
import picard.annotation.LocusFunction;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * Tags a mapped read with the genes it overlaps, the strand of those genes, the locus function (UTR/CODING/INTRONIC/...)
 * of the read on each gene, and a single summary locus function for the read.  This is the per-read work of
 * {@link TagReadWithGeneFunction}, separated from that command so it can be unit tested on its own and run on
 * multiple threads.
 *
 * <p>Thread safety: a tagger is <b>not</b> thread safe, because it accumulates strand metrics as it tags.  The
 * gene models (the {@link OverlapDetector} and the {@link Gene}s in it) are only read, so they can be shared by any
 * number of taggers.  To tag on multiple threads, give each thread its own tagger over the same detector, and combine
 * the metrics at the end with {@link TagReadWithGeneFunction.ReadTaggingMetric#merge}.
 */
public class ReadGeneFunctionTagger {

	public static final Comparator<Gene> GENE_NAME_COMPARATOR = (g1, g2) -> g1.getName().compareTo(g2.getName());

	private final OverlapDetector<Gene> geneOverlapDetector;
	private final String geneNameTag;
	private final String geneStrandTag;
	private final String geneFunctionTag;
	private final String readFunctionTag;
	private final boolean useStrandInfo;
	private final double pctRequiredOverlap;
	private final AnnotationUtils annotationUtils = AnnotationUtils.getInstance();
	private final TagReadWithGeneFunction.ReadTaggingMetric metrics = new TagReadWithGeneFunction.ReadTaggingMetric();

	/**
	 * @param geneOverlapDetector The gene models to tag reads with.
	 * @param geneNameTag Takes on the name(s) of the gene(s) the read overlaps (if any).
	 * @param geneStrandTag For each gene name, the strand of the gene.
	 * @param geneFunctionTag For each gene name, the function of the gene at the read's position.
	 * @param readFunctionTag For this read, what is the function?  There is only 1 value for this tag.
	 * @param useStrandInfo Use the read strand to select the gene(s) used for the readFunctionTag.
	 * @param pctRequiredOverlap What percent of the read must overlap a locus function for it to be included in the annotation.
	 */
	public ReadGeneFunctionTagger(final OverlapDetector<Gene> geneOverlapDetector, final String geneNameTag,
								  final String geneStrandTag, final String geneFunctionTag, final String readFunctionTag,
								  final boolean useStrandInfo, final double pctRequiredOverlap) {
		this.geneOverlapDetector = geneOverlapDetector;
		this.geneNameTag = geneNameTag;
		this.geneStrandTag = geneStrandTag;
		this.geneFunctionTag = geneFunctionTag;
		this.readFunctionTag = readFunctionTag;
		this.useStrandInfo = useStrandInfo;
		this.pctRequiredOverlap = pctRequiredOverlap;
	}

	/** The strand metrics accumulated by this tagger so far. */
	public TagReadWithGeneFunction.ReadTaggingMetric getMetrics() {
		return metrics;
	}

	/**
	 * Add functional annotations for this read.  The read is modified and returned.
	 * @param r a mapped read.
	 */
	public SAMRecord tag (final SAMRecord r) {
		Map<Gene, List<LocusFunction>> map = annotationUtils.getFunctionalDataForRead (r, geneOverlapDetector, pctRequiredOverlap);
		List<Gene> genes = new ArrayList<Gene>(map.keySet());
		// sort genes by name alphabetically to maintain some consistent ordering.
		Collections.sort(genes, GENE_NAME_COMPARATOR);


		// populate a list of genes and loci functions.  Get rid of intergenic annotations.
		List<Gene> geneList = new ArrayList<Gene>();

		// ensure the loci are sorted in the same output order as the genes!
		List<LocusFunction> functionList = new ArrayList<LocusFunction>();
		for (Gene g: genes) {
			List<LocusFunction> fa = map.get(g);
			Collections.sort(fa, Collections.reverseOrder());
			for (LocusFunction lf: fa) {
				geneList.add(g);
				functionList.add(lf);
			}
		}

		String finalGeneName = getCompoundName(geneList);
		String finalGeneStrand = getCompoundStrand(geneList);
		String finalLociFunctions=getCompoundFunctionName(functionList);

		r.setAttribute(this.geneFunctionTag, finalLociFunctions);
		r.setAttribute(this.geneNameTag, finalGeneName);
		r.setAttribute(this.geneStrandTag, finalGeneStrand);

		// filter Gene/Locus Function map down to correct strand for read tags for the summary function

		if (useStrandInfo)
			genes = getGenesConsistentWithReadStrand(genes, r);


		List<LocusFunction> finalFunctionList = new ArrayList<LocusFunction>();
		for (Gene g: genes) {
			List<LocusFunction> lf = map.get(g);
			finalFunctionList.addAll(lf);
		}
		// pick the summary locus function for genes that were consistent.
		LocusFunction f = annotationUtils.getLocusFunction(finalFunctionList, false);
		r.setAttribute(this.readFunctionTag, f.name());
		return (r);
	}

	/**
	 * If a read overlaps 1 gene, accept it.
	 * If the read overlaps two genes that are on opposite strands, then the read should be tagged with only the gene on the correct strand.
	 * If the read overlaps two genes on the same strand, then it's ambiguous and no genes should be tagged.
	 * This method also gathers metrics on the strand/gene assignment.
	 * @param genes
	 * @param r
	 * @return returns the gene the read is consistent with.
	 */
	private List<Gene> getGenesConsistentWithReadStrand(final List<Gene> genes, final SAMRecord r) {
		this.metrics.TOTAL_READS++;

		List<Gene> sameStrand = new ArrayList<Gene>();
		List<Gene> oppositeStrand = new ArrayList<Gene>();

		boolean negativeStrandRead = r.getReadNegativeStrandFlag();
		for (Gene g: genes) {
			boolean geneNegativeStrand = g.isNegativeStrand();
			// if the read and gene are on the same strand, both positive or negative...
			if ((negativeStrandRead && geneNegativeStrand) || (!negativeStrandRead && !geneNegativeStrand))
				sameStrand.add(g);
			else
				oppositeStrand.add(g);
		}

		if (sameStrand.size()==0 && oppositeStrand.size()>0) {
			this.metrics.READS_WRONG_STRAND++;
			return new ArrayList<Gene>();
		}

		// otherwise, the read is unambiguously assigned to a gene on the correct strand - the sameStrandSize must be 1 as it's not 0 and not > 1.
		if (oppositeStrand.size()>0)
			this.metrics.READ_AMBIGUOUS_GENE_FIXED++;

		this.metrics.READS_RIGHT_STRAND++;
		return sameStrand;

	}

	public static String getCompoundName (final Collection<Gene> genes) {

		if (genes.isEmpty()) return (null);


		StringBuilder result = new StringBuilder();
		Iterator<Gene> iter = genes.iterator();
		result.append(iter.next().getName());


		while (iter.hasNext()) {
			result.append(",");
			result.append(iter.next().getName());
		}

		return (result.toString());
	}

	public static String getCompoundFunctionName (final Collection<LocusFunction> functions) {

		if (functions.isEmpty()) return (null);

		StringBuilder result = new StringBuilder();
		Iterator<LocusFunction> iter = functions.iterator();
		result.append(iter.next().name());


		while (iter.hasNext()) {
			result.append(",");
			result.append(iter.next().name());
		}

		return (result.toString());
	}


	public static String getCompoundStrand (final Collection<Gene> genes) {

		if (genes.isEmpty()) return (null);

		StringBuilder result = new StringBuilder();
		Iterator<Gene> iter = genes.iterator();
		result.append(Utils.negativeStrandToString(iter.next().isNegativeStrand()));

		while (iter.hasNext()) {
			result.append(",");
			result.append(Utils.negativeStrandToString(iter.next().isNegativeStrand()));
		}

		return (result.toString());
	}
}
