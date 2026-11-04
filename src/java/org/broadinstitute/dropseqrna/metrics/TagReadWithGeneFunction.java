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
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.dropseqrna.annotation.GeneAnnotationReader;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.BatchingIterable;
import org.broadinstitute.dropseqrna.utils.OrderedConcurrentMapper;
import org.broadinstitute.dropseqrna.utils.SamHeaderUtil;
import picard.annotation.Gene;
import picard.annotation.LocusFunction;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.nio.PicardHtsPath;

import java.io.File;
import java.util.*;

@CommandLineProgramProperties(
        summary = "A special case tagger.  Tags reads that are exonic for the gene name of the overlapping exon.  This is done specifically to solve the case where a read" +
			"may be tagged with a gene and an exon, but the read may not be exonic for all genes tagged.  This limits the list of genes to only those where the read overlaps the exon and the gene." +
			"Reads that overlap multiple genes are assigned to the gene that shares the strand with the read.  If that assignment is ambiguous (2 or more genes share the strand of the read), then the read is not assigned any genes.",
        oneLineSummary = "Tags gene/exons in a strand-specific way, adds locus function type. Used before running digital expression.",
        programGroup = DropSeq.class
)
public class TagReadWithGeneFunction extends CommandLineProgram {

	private final Log log = Log.getInstance(TagReadWithGeneFunction.class);
	private ProgressLogger pl = new ProgressLogger(log);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze")
	public PicardHtsPath INPUT;

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

	@Argument(doc="What fraction of the read must overlap a locus function to be included in the annotation.  If this parameter is set to 0, then if any" +
			"bases are overlap an annotation it is included in the output.  StarSolo / CellRanger would set this to 50.0, then only" +
			"functional annotations with > 50% overlap would be included.  This forces a gene to have at most one annotation (coding, intronic, etc) instead" +
			"of multiple annotations (coding + intronic). Set to 0 to reproduce the DropSeq standard, or 50 to emulate StarSolo.")
	public double PCT_REQUIRED_OVERLAP=0;

	@Argument(doc="The number of threads used to tag reads.  Reads are read and written in their original order, and the output is" +
			" identical regardless of this value.  With 1 (the default) reads are tagged one at a time on the main thread.  With more than 1, " +
			"that many threads tag reads in parallel, in addition to the thread that reads the input and the thread that writes the output.  " +
			"More threads help most when the input and output are uncompressed, since otherwise reading and writing the BAM is the bottleneck.")
	public int NUM_THREADS=1;

	@Hidden
	@Argument(doc="The number of reads handed to a tagging thread at a time when NUM_THREADS > 1.  The default is fine for almost all data.")
	public int BATCH_SIZE=1000;

	// @Option(doc="Allow a read to span the exons of multiple genes.  If set to true, the gene name will be set to all of the gene/exons the read spans.  In that case, the gene names will be comma separated.")
	private boolean ALLOW_MULTI_GENE_READS=false;

	// in-flight batches allowed per tagging thread when NUM_THREADS > 1.
	private static final int MAX_BATCHES_PER_THREAD = 4;

	private ReadTaggingMetric metrics = new ReadTaggingMetric();

	public TagReadWithGeneFunction() {
		VALIDATION_STRINGENCY = ValidationStringency.LENIENT;
	}

	@Override
	protected String[] customCommandLineValidation() {
		final List<String> errors = new ArrayList<>();
		if (NUM_THREADS < 1) errors.add("NUM_THREADS must be at least 1: " + NUM_THREADS);
		if (BATCH_SIZE < 1) errors.add("BATCH_SIZE must be at least 1: " + BATCH_SIZE);
		return errors.isEmpty() ? super.customCommandLineValidation() : errors.toArray(new String[0]);
	}

	@Override
	protected int doWork() {
		IOUtil.assertFileIsReadable(this.INPUT.toPath());
		IOUtil.assertFileIsReadable(this.ANNOTATIONS_FILE);
		if (this.SUMMARY!=null) IOUtil.assertFileIsWritable(this.SUMMARY);
		IOUtil.assertFileIsWritable(this.OUTPUT);

		final SamReader inputSam = SamReaderFactory.makeDefault().open(INPUT.toPath());

		SAMFileHeader header = inputSam.getFileHeader();
		SamHeaderUtil.addPgRecord(header, this);
		SAMSequenceDictionary bamDict = header.getSequenceDictionary();

        final OverlapDetector<Gene> geneOverlapDetector = GeneAnnotationReader.loadAnnotationsFile(ANNOTATIONS_FILE, bamDict,
				VALIDATION_STRINGENCY);
        SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT);

        if (NUM_THREADS == 1)
        	tagReads(inputSam, writer, geneOverlapDetector);
        else
        	tagReadsConcurrently(inputSam, writer, geneOverlapDetector);

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
	private ReadGeneFunctionTagger newTagger(final OverlapDetector<Gene> geneOverlapDetector) {
		return new ReadGeneFunctionTagger(geneOverlapDetector, GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG,
				READ_FUNCTION_TAG, USE_STRAND_INFO, PCT_REQUIRED_OVERLAP);
	}

	/**
	 * Tag the reads one at a time on this thread, writing them in the order they were read.
	 */
	private void tagReads(final SamReader inputSam, final SAMFileWriter writer, final OverlapDetector<Gene> geneOverlapDetector) {
		final ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		for (SAMRecord r: inputSam) {
			pl.record(r);
			if (!r.getReadUnmappedFlag())
				tagger.tag(r);
			writer.addAlignment(r);
		}
		this.metrics.merge(tagger.getMetrics());
	}

	/**
	 * Tag the reads on NUM_THREADS threads, writing them in the order they were read.
	 *
	 * The calling thread reads the input and hands batches of reads to the tagging threads.  Each tagging thread has its
	 * own tagger (which accumulates its own metrics), over the same read-only gene models.  A single thread writes
	 * the tagged batches in the original order.  OrderedConcurrentMapper counts that writing thread as one of its
	 * threads, hence NUM_THREADS + 1.
	 */
	private void tagReadsConcurrently(final SamReader inputSam, final SAMFileWriter writer, final OverlapDetector<Gene> geneOverlapDetector) {
		final List<ReadGeneFunctionTagger> taggers = Collections.synchronizedList(new ArrayList<>());
		final Iterator<SAMRecord> progressLoggingIterator = new Iterator<>() {
			private final Iterator<SAMRecord> it = inputSam.iterator();
			@Override
			public boolean hasNext() { return it.hasNext(); }
			@Override
			public SAMRecord next() {
				final SAMRecord r = it.next();
				pl.record(r);
				return r;
			}
		};

		new OrderedConcurrentMapper.Builder<List<SAMRecord>, ReadGeneFunctionTagger, List<SAMRecord>>()
				.withIterable(new BatchingIterable<>(progressLoggingIterator, BATCH_SIZE))
				.withNumThreads(NUM_THREADS + 1)
				// bounds the batches that have been read but not yet written.  OrderedConcurrentMapper's default
				// counts items, and is meant for single records, not batches of them.
				.withMaxItems(MAX_BATCHES_PER_THREAD * NUM_THREADS)
				.withPerThreadFactory(() -> {
					final ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
					taggers.add(tagger);
					return tagger;
				})
				.withMapper((final List<SAMRecord> batch, final ReadGeneFunctionTagger tagger) -> {
					for (final SAMRecord r : batch)
						if (!r.getReadUnmappedFlag())
							tagger.tag(r);
					return batch;
				})
				.withConsumer(batch -> batch.forEach(writer::addAlignment))
				.build()
				.execute();

		// every tagging thread has finished, so the taggers' metrics are safe to read.
		for (final ReadGeneFunctionTagger tagger : taggers)
			this.metrics.merge(tagger.getMetrics());
	}

	/**
	 * Add functional annotations for this read.  Convenience for tagging a single read.
	 * @param r a mapped read
	 * @param geneOverlapDetector the gene models
	 * @param allowMultiReadGenes unused
	 * @return The read.
	 */
	public SAMRecord setAnnotations (final SAMRecord r, final OverlapDetector<Gene> geneOverlapDetector, final boolean allowMultiReadGenes) {
		final ReadGeneFunctionTagger tagger = newTagger(geneOverlapDetector);
		tagger.tag(r);
		this.metrics.merge(tagger.getMetrics());
		return r;
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

	public static final Comparator<Gene> GENE_NAME_COMPARATOR = ReadGeneFunctionTagger.GENE_NAME_COMPARATOR;


    public static final Comparator<LocusFunction> LOCUS_FUNCTION_COMPARATOR =  new Comparator<LocusFunction>() {
        @Override
		public int compare(final LocusFunction g1, final LocusFunction g2) {
            return (g1.compareTo(g2));
        }
    };

	public static String getCompoundName (final Collection<Gene> genes) {
		return ReadGeneFunctionTagger.getCompoundName(genes);
	}

	public static String getCompoundFunctionName (final Collection<LocusFunction> functions) {
		return ReadGeneFunctionTagger.getCompoundFunctionName(functions);
	}

	public static String getCompoundStrand (final Collection<Gene> genes) {
		return ReadGeneFunctionTagger.getCompoundStrand(genes);
	}

	public static class ReadTaggingMetric extends MetricBase {
		public int TOTAL_READS=0;
		public int READS_WRONG_STRAND=0;
		public int READS_RIGHT_STRAND=0;
		public int READ_AMBIGUOUS_GENE_FIXED=0;
		public int AMBIGUOUS_READS_REJECTED=0;

		/** Add the counts of another metric to this one. */
		public void merge(final ReadTaggingMetric other) {
			this.TOTAL_READS += other.TOTAL_READS;
			this.READS_WRONG_STRAND += other.READS_WRONG_STRAND;
			this.READS_RIGHT_STRAND += other.READS_RIGHT_STRAND;
			this.READ_AMBIGUOUS_GENE_FIXED += other.READ_AMBIGUOUS_GENE_FIXED;
			this.AMBIGUOUS_READS_REJECTED += other.AMBIGUOUS_READS_REJECTED;
		}

		@Override
		public String toString () {
			return("TOTAL READS [" + this.TOTAL_READS + "] CORRECT_STRAND [" + this.READS_RIGHT_STRAND +"]  WRONG_STRAND [" + this.READS_WRONG_STRAND +"] AMBIGUOUS_STRAND_FIXED [" + this.READ_AMBIGUOUS_GENE_FIXED +"] AMBIGUOUS REJECTED READS ["+ this.AMBIGUOUS_READS_REJECTED+"]" );
		}
	}

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new TagReadWithGeneFunction().instanceMain(args));
	}

}
