/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.dropseqrna.utils.alignmentcomparison;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.stream.Collectors;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.StringInterner;
import org.broadinstitute.dropseqrna.utils.alignmentcomparison.QueryNameJointIterator.JointResult;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.MapQualityFilteredIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamRecordSortingIteratorFactory;

import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.StringUtil;
import picard.cmdline.CommandLineProgram;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.Argument;

/**
 *
 * @author nemesh
 *
 */
@CommandLineProgramProperties(summary = "Compare two dropseq alignments across different sets of annotations or genome builds to see what genes present in the original set are lost in the second set.  "
		+ "The input data should be the same set of data aligned to different references, on either the same genome build or different builds.  If the second input is aligned with secondary reads emitted during alignment,"
		+ "all contigs for a mapping for non-unique contigs will be listed.  This can be helpful to detect if additional contigs [such as alternative haplotypes] create mapping problems that make some subset of genes invisible to the dropseq toolkit.",
oneLineSummary = "Compare two alignments",
programGroup = DropSeq.class)
public class CompareDropSeqAlignments extends CommandLineProgram {

	private static final Log log = Log.getInstance(CompareDropSeqAlignments.class);

	@Argument(doc = "The input SAM or BAM file to analyze.  If coordinate sorted this will save time, but is not required.")
	public File INPUT_1;

	@Argument(doc = "The comparison input SAM or BAM file to analyze.  If coordinate sorted this will save time, but is not required.")
	public File INPUT_2;

	@Argument(doc="Output file that maps the contig the read uniquely mapped to in INPUT_1, and the contig the read mapped to in INPUT_2, with reads partitioned into groups that did/did not remained uniquely mapped.  This supports zipped formats like gz and bz2.", optional=true)
	public File CONTIG_REPORT;

	@Argument(doc="Output file that maps the gene the read uniquely mapped to in INPUT_1, and the gene the read mapped to in INPUT_2, with reads partitioned into groups that did/did not remained uniquely mapped.  This supports zipped formats like gz and bz2.", optional=true)
	public File GENE_REPORT;

	@Argument(doc="Tag to extract", optional=true)
	public String GENE_EXON_TAG="GE";

	@Argument(doc="The map quality for a read to be considered uniquely mapped.")
	public int READ_QUALITY=10;

	@Argument (doc="Trim this string from the contig names of both BAMs to make contig names comparable.  This is useful when one alignment strategy calls the first contig 'chr1' and the second strategy '1'")
	public String TRIM_CONTIG_STRING="chr";

	private StringInterner stringInterner = new StringInterner();
	private final String noGeneTag="NO_GENE";

	// groups of contigs: alt, random, HLA, decoy
	@Override
	protected int doWork() {
		ObjectCounter<ContigResult> contigResults = null;
		// key = original gene name, value = collection of mappings.
		Map<String, GeneResult> geneResults = null;
		IOUtil.assertFileIsReadable(INPUT_1);
		IOUtil.assertFileIsReadable(INPUT_2);
		if (this.CONTIG_REPORT==null && this.GENE_REPORT==null) log.error ("What's the point of running me if there's no output defined?");
		if (this.CONTIG_REPORT!=null) {
			IOUtil.assertFileIsWritable(this.CONTIG_REPORT);
			contigResults = new ObjectCounter<>();
		}
		if (this.GENE_REPORT!=null) {
			IOUtil.assertFileIsWritable(this.GENE_REPORT);
			if (this.GENE_EXON_TAG==null) {
				log.error("If GENE_REPORT is specified, a GENE_EXON_TAG must also be specified");
				return 1;
			}

			geneResults = new HashMap<>();
		}

		PeekableIterator<List<SAMRecord>> oIter = getReadIterator(this.INPUT_1, this.READ_QUALITY);
		PeekableIterator<List<SAMRecord>> nIter = getReadIterator(this.INPUT_2, null);

		QueryNameJointIterator qnji = new QueryNameJointIterator(oIter, nIter);

		int counter =0;
		while (qnji.hasNext()) {
			if (counter%1000000==0) log.info("Number of reads processed [" + counter + "]");
			JointResult jr = qnji.next();
			counter++;

			List<SAMRecord> r1 = jr.getOne();
			List<SAMRecord> r2 = jr.getTwo();

			if (contigResults!=null) {
				ContigResult cr = evaluateByContig(r1, r2);
				contigResults.increment(cr);
			}

			if (geneResults!=null)
				geneResults=evaluateByGene(r1, r2, geneResults, this.GENE_EXON_TAG);

		}

		if (this.CONTIG_REPORT!=null) writeContigReport(this.CONTIG_REPORT, contigResults);
		if (this.GENE_REPORT!=null) writeGeneReport (this.GENE_REPORT, geneResults);
		return 0;
	}


	/**
	 * Compares the geneExonTag of the first read (which should be a uniquely mapped read, but may not have a gene-exon tag) to the
	 * second read(s).
	 * @param r1
	 * @param r2
	 * @param geneResults
	 * @param geneExonTag
	 * @return
	 */
	private Map<String, GeneResult> evaluateByGene (final List<SAMRecord> r1, final List<SAMRecord> r2, final Map<String, GeneResult> geneResults, final String geneExonTag) {
		if (!validateReadSetSize(r1, r2)) return geneResults;
		SAMRecord rr1 = r1.get(0);
		String originalGene = rr1.getStringAttribute(geneExonTag);
		// if there's no original gene mapping, skip.
		if (originalGene==null) return geneResults;
		// otherwise fetch the original gene if it exists, or add it to the map.
		originalGene=this.stringInterner.intern(originalGene);
		GeneResult gr = geneResults.get(originalGene);
		// log.info(r1.get(0).getReadName());
		if (gr==null)
			gr=new GeneResult(originalGene, rr1.getContig(), this.noGeneTag);

		// need to replace NULL with "NO_GENE" instead of filtering result.
		Function<String,String> fixGeneName = (final String geneName) -> {
			if (geneName!=null) return geneName;
			else return this.noGeneTag;
		};

		Collection<String> genesNew = r2.stream().map(x-> x.getStringAttribute(geneExonTag)).map(fixGeneName).map(x-> stringInterner.intern(x)).collect(Collectors.toList());
		// make this collection unique genes...
		List<String> geneNewUnique = new ArrayList<>(new TreeSet<>(genesNew));

		// get contigs.
		Collection<String> contigsNew = r2.stream().map(x-> x.getContig()).filter(x->x!=null).map(x -> x.replaceAll(this.TRIM_CONTIG_STRING, "")).map(x-> stringInterner.intern(x)).collect(Collectors.toList());
		// get unique contigs.
		List<String> contigsNewUnique = new ArrayList<>(new TreeSet<>(contigsNew));
		gr.addMapping(geneNewUnique, contigsNewUnique, r2.size());

		geneResults.put(originalGene, gr);
		return geneResults;
	}







	private ContigResult evaluateByContig (final List<SAMRecord> r1, final List<SAMRecord> r2) {
		if (!validateReadSetSize(r1, r2)) return null;
		String contigOne = stringInterner.intern(r1.get(0).getContig().replaceAll(this.TRIM_CONTIG_STRING, ""));
		Collection<String> contigsNew = r2.stream().map(x-> x.getContig()).map(x -> x.replaceAll(this.TRIM_CONTIG_STRING, "")).map(x-> stringInterner.intern(x)).collect(Collectors.toList());
		List<String> contigsNewUnique = new ArrayList<>(new TreeSet<>(contigsNew));
		ContigResult r = new ContigResult(contigOne, contigsNewUnique, contigsNewUnique.size()==1);
		return (r);
	}

	private PeekableIterator<List<SAMRecord>> getReadIterator (final File bamFile, final Integer readQuality) {
        Iterator<SAMRecord> iter = getQueryNameSortedData(bamFile);
		// filter out unmapped reads.
		iter = new UnmappedReadFilter(iter);
		// optionally, filter out reads below a map quality threshold.
		if (readQuality!=null) iter = new MapQualityFilteredIterator(iter, readQuality, false).iterator();
		final GroupingIterator<SAMRecord> groupingIterator = new GroupingIterator<>(iter, READ_NAME_COMPARATOR);
		PeekableIterator<List<SAMRecord>> peekable = new PeekableIterator<>(groupingIterator);
		return peekable;
	}

	private Iterator<SAMRecord> getQueryNameSortedData (final File bamFile) {
		SamReader reader = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE).open(bamFile);
		if (reader.getFileHeader().getSortOrder().equals(SortOrder.queryname))
			return reader.iterator();
		log.info("Input SAM/BAM not in queryname order, sorting...");
        final ProgressLogger progressLogger = new ProgressLogger(log, 1000000, "Sorting reads in query name order");
        final CloseableIterator<SAMRecord> result = SamRecordSortingIteratorFactory.create(reader.getFileHeader(), reader.iterator(), READ_NAME_COMPARATOR, progressLogger);
        log.info("Sorting finished.");
        return result;
	}

	private void writeGeneReport (final File outFile, final Map<String, GeneResult> geneResults) {
		PrintStream writer = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(outFile));
		List<String> header = new ArrayList<>();
		header.add("INPUT_1="+this.INPUT_1.getAbsolutePath());
		header.add("INPUT_2="+this.INPUT_2.getAbsolutePath());
		header.add("READ_QUALITY="+this.READ_QUALITY);
		header.add("GENE_EXON_TAG="+this.GENE_EXON_TAG);
		String h = StringUtils.join(header, "\t");
		writer.print("#");
		writer.println(h);

		String [] colNames = {"GENE", "CONTIG", "ORIGINAL_READS", "SAME_MAPPING",
				"NO_GENE", "SAME_GENE_MULTIREAD",
				"DIFFERENT_GENE", "DIFFERENT_GENE_MULTIREAD","MULTI_GENE_MAPPING",
				"OTHER_UNIQUE_GENES", "OTHER_MULTIMAPPING_GNES", "OTHER_CONTIGS"};
		writer.println(StringUtil.join("\t", colNames));


		List<String> keys = new ArrayList<>(geneResults.keySet());
		Collections.sort(keys);

		for (String key: keys) {
			GeneResult gr = geneResults.get(key);
			String [] line = {gr.getOriginalGene(),  gr.getOriginalContig(), Integer.toString(gr.getCountOriginalReads()), Integer.toString(gr.getCountSameMapping()),
					Integer.toString(gr.getCountIntronicOrIntergenic()), Integer.toString(gr.getCountSameGeneMapsNonUniqueCount()),
					Integer.toString(gr.getCountDifferentUniqueGene()), Integer.toString(gr.getCountDifferentGeneNonUniqueCount()), Integer.toString(gr.getCountMultiGeneMappingCount()),
					StringUtil.join(",", gr.getUniqueMapOtherGene()), StringUtil.join(",", gr.getNonUniqueMapOtherGene()), StringUtil.join(",", gr.getOtherContigs())};
			writer.println(StringUtil.join("\t", line));
		}
		writer.close();
	}

	private void writeContigReport (final File outFile, final ObjectCounter<ContigResult> contigResults) {
		PrintStream writer = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(outFile));
		List<String> header = new ArrayList<>();
		header.add("INPUT_1="+this.INPUT_1.getAbsolutePath());
		header.add("INPUT_2="+this.INPUT_2.getAbsolutePath());
		header.add("READ_QUALITY="+this.READ_QUALITY);
		header.add("TRIM_CONTIG_STRING="+this.TRIM_CONTIG_STRING);
		String h = StringUtils.join(header, "\t");
		writer.print("#");
		writer.println(h);

		String [] colNames = {"ORIGINAL_CONTIG", "NEW_CONTIGS", "NUM_READS", "MAPPED_UNIQUELY"};
		writer.println(StringUtil.join("\t", colNames));
		List<ContigResult> crList = contigResults.getKeysOrderedByCount(true);
		for (ContigResult cr: crList) {
			String [] body = {cr.getOldContig(), StringUtil.join(",", cr.getNewContigs()), Integer.toString(contigResults.getCountForKey(cr)), Boolean.toString(cr.isNewReadMapsUniquely())};
			writer.println(StringUtil.join("\t", body));
		}
		writer.close();
	}


	private boolean validateReadSetSize (final List<SAMRecord> r1, final List<SAMRecord> r2) {
		// handle multi-mapping inputs by skipping.
		if (r1.size()>1) return false;
		if (r1.size()==0) throw new IllegalStateException ("First read set size = 0.  This should never happen");
		if (r2.size()==0) throw new IllegalStateException ("Second read set size = 0.  This should never happen");
		return true;
	}

	static final Comparator<SAMRecord> READ_NAME_COMPARATOR =  new Comparator<SAMRecord>() {
        private final SAMRecordQueryNameComparator comp = new SAMRecordQueryNameComparator();
        @Override
		public int compare(final SAMRecord s1, final SAMRecord s2) {
            return comp.fileOrderCompare(s1, s2);
        }
    };

    /**
     * Filters out unmapped reads.
     * @author nemesh
     *
     */
    public class UnmappedReadFilter extends FilteredIterator<SAMRecord> {
    	public UnmappedReadFilter(final Iterator<SAMRecord> underlyingIterator) {
            super(underlyingIterator);
    	}
        @Override
        public boolean filterOut(final SAMRecord r) {
        	return r.getReadUnmappedFlag();
        }
    }


	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new CompareDropSeqAlignments().instanceMain(args));
	}
}


