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
package org.broadinstitute.dropseqrna.metrics;

import htsjdk.samtools.*;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.annotation.AnnotationUtils;
import org.broadinstitute.dropseqrna.annotation.GeneAnnotationReader;
import org.broadinstitute.dropseqrna.barnyard.Utils;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.SamHeaderUtil;
import picard.annotation.Gene;
import picard.annotation.LocusFunction;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.util.*;

@CommandLineProgramProperties(
        summary = "A special case tagger.  Tags reads that are exonic for the gene name of the overlapping exon.  This is done specifically to solve the case where a read" +
			"may be tagged with a gene and an exon, but the read may not be exonic for all genes tagged.  This limits the list of genes to only those where the read overlaps the exon and the gene." +
			"Reads that overlap multiple genes are assigned to the gene that shares the strand with the read.  If that assignment is ambiguous (2 or more genes share the strand of the read), then the read is not assigned any genes.",
        oneLineSummary = "Tags gene/exons in a strand-specific way, adds locus function type. Used before running digital expression.",
        programGroup = DropSeq.class
)
public class TagReadWithGeneExon extends CommandLineProgram {

	private final Log log = Log.getInstance(TagReadWithGeneExon.class);
	private ProgressLogger pl = new ProgressLogger(log);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze")
	public File INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output BAM, written with new Gene/Exon tag")
	public File OUTPUT;

	@Argument(doc = "The strand specific summary info", optional=true)
	public File SUMMARY=null;

	/*
	@Option(doc = "The tag name to use.  When there are multiple genes, they will be comma seperated.")
	public String TAG="GE";

	@Argument(doc="The strand of the gene(s) the read overlaps.  When there are multiple genes, they will be comma seperated.")
	public String STRAND_TAG="GS";

	@Argument(doc = "The functional annotation for the read")
	public String FUNCTION_TAG="XF";
	*/

	@Argument(doc="Gene Name tag.  Takes on the gene name this read overlaps (if any)")
	public String GENE_NAME_TAG="gn";

	@Argument(doc="Gene Strand tag.  For a given gene name <GENE_NAME_TAG>, this is the strand of the gene.")
	public String GENE_STRAND_TAG="gs";

	@Argument(doc="Gene Function tag.  For a given gene name <GENE_NAME_TAG>, this is the function of the gene at this read's position: UTR/CODING/INTRONIC/...")
	public String GENE_FUNCTION_TAG="gf";

	@Argument(doc="READ functional annotation tag.  For this read, what is the function?  This only looks at reads on the right strand, and prioritizes coding over intron over intergenic.  There is only 1 value for this tag.")
	public String READ_FUNCTION_TAG="XF";

	@Argument(doc="The annotations set to use to label the read.  This can be a GTF or a refFlat file.")
	public File ANNOTATIONS_FILE;

	@Argument(doc="Use strand info to determine what gene to assign the read to.  If this is on, reads can be assigned to a maximum one one gene.  This is used for the READ_FUNCTION_TAG output only.")
	public boolean USE_STRAND_INFO=true;

	// @Option(doc="Allow a read to span the exons of multiple genes.  If set to true, the gene name will be set to all of the gene/exons the read spans.  In that case, the gene names will be comma separated.")
	private boolean ALLOW_MULTI_GENE_READS=false;

	private ReadTaggingMetric metrics = new ReadTaggingMetric();

	@Override
	protected int doWork() {
		IOUtil.assertFileIsReadable(this.INPUT);
		IOUtil.assertFileIsReadable(this.ANNOTATIONS_FILE);
		if (this.SUMMARY!=null) IOUtil.assertFileIsWritable(this.SUMMARY);
		IOUtil.assertFileIsWritable(this.OUTPUT);

		SamReader inputSam = SamReaderFactory.makeDefault().open(INPUT);

		SAMFileHeader header = inputSam.getFileHeader();
		SamHeaderUtil.addPgRecord(header, this);
		SAMSequenceDictionary bamDict = header.getSequenceDictionary();

        final OverlapDetector<Gene> geneOverlapDetector = GeneAnnotationReader.loadAnnotationsFile(ANNOTATIONS_FILE, bamDict);
        SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT);

        for (SAMRecord r: inputSam) {
        	pl.record(r);

        	if (!r.getReadUnmappedFlag())
				// r=	setGeneExons(r, geneOverlapDetector, this.ALLOW_MULTI_GENE_READS);
        		r= setAnnotations(r, geneOverlapDetector, this.ALLOW_MULTI_GENE_READS);
        	writer.addAlignment(r);
        }

		CloserUtil.close(inputSam);
		writer.close();
		// if (this.USE_STRAND_INFO) log.info(this.metrics.toString());
		if (SUMMARY==null) return 0;

		//process summary
		MetricsFile<ReadTaggingMetric, Integer> outFile = new MetricsFile<ReadTaggingMetric, Integer>();
		outFile.addMetric(this.metrics);
		outFile.write(this.SUMMARY);
		return 0;

	}

	/*
	public SAMRecord setGeneExons (final SAMRecord r, final OverlapDetector<Gene> geneOverlapDetector, final boolean allowMultiGeneReads) {
		Map<Gene, LocusFunction> map = AnnotationUtils.getInstance().getLocusFunctionForReadByGene(r, geneOverlapDetector);
		Set<Gene> exonsForRead = AnnotationUtils.getInstance().getConsistentExons (r, map.keySet(), allowMultiGeneReads);

		List<Gene> genes = new ArrayList<Gene>();

 		for (Gene g: exonsForRead) {

			LocusFunction f = map.get(g);
			if (f==LocusFunction.CODING || f==LocusFunction.UTR)
				genes.add(g);
		}
		LocusFunction f = AnnotationUtils.getInstance().getLocusFunction(map.values());

		if (USE_STRAND_INFO)
			genes = getGenesConsistentWithReadStrand(genes, r);

		if (genes.size()>1 && this.ALLOW_MULTI_GENE_READS==false)
			log.error("There should only be 1 gene assigned to a read for DGE purposes.");

		String finalGeneName = getCompoundName(genes);
		String finalGeneStrand = getCompoundStrand(genes);

		if (f!=null)
			r.setAttribute(this.FUNCTION_TAG, f.toString());
		if (finalGeneName!=null && finalGeneStrand!=null) {
			r.setAttribute(this.TAG, finalGeneName);
			r.setAttribute(this.STRAND_TAG, finalGeneStrand);
		} else {
			r.setAttribute(this.TAG, null);
			r.setAttribute(this.STRAND_TAG, null);
		}
		return (r);
	}
	*/
	/**
	 * Add functional annotations for this read.
	 * @param r
	 * @param geneOverlapDetector
	 * @return
	 */
	public SAMRecord setAnnotations (final SAMRecord r, final OverlapDetector<Gene> geneOverlapDetector, final boolean allowMultiReadGenes) {
		Map<Gene, List<LocusFunction>> map = AnnotationUtils.getInstance().getFunctionalDataForRead (r, geneOverlapDetector);
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

		r.setAttribute(this.GENE_FUNCTION_TAG, finalLociFunctions);
		r.setAttribute(this.GENE_NAME_TAG, finalGeneName);
		r.setAttribute(this.GENE_STRAND_TAG, finalGeneStrand);

		// filter Gene/Locus Function map down to correct strand for read tags for the summary function

		if (USE_STRAND_INFO)
			genes = getGenesConsistentWithReadStrand(genes, r);


		List<LocusFunction> finalFunctionList = new ArrayList<LocusFunction>();
		for (Gene g: genes) {
			List<LocusFunction> lf = map.get(g);
			finalFunctionList.addAll(lf);
		}
		// pick the summary locus function for genes that were consistent.
		LocusFunction f = AnnotationUtils.getInstance().getLocusFunction(finalFunctionList, false);
		r.setAttribute(this.READ_FUNCTION_TAG, f.name());
		return (r);
	}

	/**
	 * Get a list of genes that are not intergenic.  Don't record intergenic annotations, as those are genes that fall inbetween alignment blocks.
	 * @param map
	 * @return
	 */
	private Set<Gene> getNonIntergenicGenes (final Map<Gene, LocusFunction> map) {
		Set<Gene> result = new HashSet<Gene>();
		for (Gene g: map.keySet()) {
			LocusFunction lf = map.get(g);
			if (lf!=LocusFunction.INTERGENIC)
				result.add(g);
		}
		return (result);
	}

	public static final Comparator<Gene> GENE_NAME_COMPARATOR =  new Comparator<Gene>() {
        @Override
		public int compare(final Gene g1, final Gene g2) {
            return (g1.getName().compareTo(g2.getName()));
        }
    };

    public static final Comparator<LocusFunction> LOCUS_FUNCTION_COMPARATOR =  new Comparator<LocusFunction>() {
        @Override
		public int compare(final LocusFunction g1, final LocusFunction g2) {
            return (g1.compareTo(g2));
        }
    };

    // LocusFunction

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

		/**
		if (sameStrand.size()>1) {
			this.metrics.AMBIGUOUS_READS_REJECTED++;
			return new ArrayList<Gene>();
		}
		*/
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
		result.append(Utils.strandToString(iter.next().isPositiveStrand()));

		while (iter.hasNext()) {
			result.append(",");
			result.append(Utils.strandToString(iter.next().isPositiveStrand()));
		}

		return (result.toString());
	}



	public class ReadTaggingMetric extends MetricBase {
		public int TOTAL_READS=0;
		public int READS_WRONG_STRAND=0;
		public int READS_RIGHT_STRAND=0;
		public int READ_AMBIGUOUS_GENE_FIXED=0;
		public int AMBIGUOUS_READS_REJECTED=0;

		@Override
		public String toString () {
			return("TOTAL READS [" + this.TOTAL_READS + "] CORRECT_STRAND [" + this.READS_RIGHT_STRAND +"]  WRONG_STRAND [" + this.READS_WRONG_STRAND +"] AMBIGUOUS_STRAND_FIXED [" + this.READ_AMBIGUOUS_GENE_FIXED +"] AMBIGUOUS REJECTED READS ["+ this.AMBIGUOUS_READS_REJECTED+"]" );
		}
	}

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new TagReadWithGeneExon().instanceMain(args));
	}

}
