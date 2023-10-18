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
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.StringInterner;
import org.broadinstitute.dropseqrna.utils.alignmentcomparison.QueryNameJointIterator.JointResult;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.MapQualityFilteredIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamRecordSortingIteratorFactory;

import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.StringUtil;
import picard.cmdline.CommandLineProgram;

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
@groovy.transform.Generated
public class CompareDropSeqAlignments  extends CommandLineProgram {

	private static final Log log = Log.getInstance(CompareDropSeqAlignments.class);

	@Argument(doc = "The input SAM or BAM file to analyze.  If query name sorted this will save time, but is not required.", minElements = 1)
	public List<File> INPUT_1;

	@Argument(doc = "The comparison input SAM or BAM file to analyze.  If query name sorted this will save time, but is not required.", minElements = 1)
	public List<File> INPUT_2;

	@Argument(doc="Output file that maps the contig the read uniquely mapped to in INPUT_1, and the contig the read mapped to in INPUT_2, with reads partitioned into groups that did/did not remained uniquely mapped.  This supports zipped formats like gz and bz2.", optional=true)
	public File CONTIG_REPORT;
	
	@Argument(doc="Output each read that disagrees on mapping by contig/position.  Emits the read name, location, and map quality of the read in both alignments.  Also emits the TAG_NAME and values for those tags for the read in each alignment.", optional=true)
	public File READ_REPORT;

	@Argument(doc="A list of 1 or more tags who's values can be reported for read in the READ_REPORT.", optional=true)
	public List<String> TAG_NAMES;
	
	// @Argument(doc="Output file that maps the gene the read uniquely mapped to in INPUT_1, and the gene the read mapped to in INPUT_2, with reads partitioned into groups that did/did not remained uniquely mapped.  This supports zipped formats like gz and bz2.", optional=true)
	// public File GENE_REPORT;

	// @Argument(doc="Tag to extract", optional=true)
	// public String GENE_EXON_TAG="GE";

	@Argument(doc="The map quality for a read to be considered uniquely mapped.")
	public int READ_QUALITY=10;

	@Argument (doc="Trim this string from the contig names of both BAMs to make contig names comparable.  This is useful when one alignment strategy calls the first contig 'chr1' and the second strategy '1'")
	public String TRIM_CONTIG_STRING="";

	private StringInterner stringInterner = new StringInterner();
	final static String noGeneTag="NO_GENE";

	// groups of contigs: alt, random, HLA, decoy
	@Override
	protected int doWork() {
		ObjectCounter<ContigResult> contigResults = null;
		// key = original gene name, value = collection of mappings.
		// Map<String, GeneResult> geneResults = null;
		INPUT_1=FileListParsingUtils.expandFileList(INPUT_1);
		INPUT_2=FileListParsingUtils.expandFileList(INPUT_2);
		INPUT_1.stream().forEach(x-> IOUtil.assertFileIsReadable(x));
		INPUT_2.stream().forEach(x-> IOUtil.assertFileIsReadable(x));
		// if (this.CONTIG_REPORT==null && this.GENE_REPORT==null) log.error ("What's the point of running me if there's no output defined?");
		if (this.CONTIG_REPORT!=null) {
			IOUtil.assertFileIsWritable(this.CONTIG_REPORT);
			contigResults = new ObjectCounter<>();
		}
		
		PrintStream readWriter=null;
		if (this.READ_REPORT!=null) {
			IOUtil.assertFileIsWritable(this.READ_REPORT);
			readWriter = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(this.READ_REPORT));	
			writeReadHeader(this.TAG_NAMES, readWriter);
		}
		
		/*
		if (this.GENE_REPORT!=null) {
			IOUtil.assertFileIsWritable(this.GENE_REPORT);
			if (this.GENE_EXON_TAG==null) {
				log.error("If GENE_REPORT is specified, a GENE_EXON_TAG must also be specified");
				return 1;
			}

			geneResults = new HashMap<>();
		}
		*/
		
		log.info("Reading INPUT_1");
		PeekableIterator<List<SAMRecord>> oIter = getReadIterator(this.INPUT_1, this.READ_QUALITY, true);
		log.info("Reading INPUT_2");
		PeekableIterator<List<SAMRecord>> nIter = getReadIterator(this.INPUT_2, null, false);

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
			if (readWriter!=null) {
				writeReadReportLine(r1, r2, this.TAG_NAMES, readWriter);
			}
			/*
			if (geneResults!=null)
				geneResults=evaluateByGene(r1, r2, geneResults, this.GENE_EXON_TAG);
			*/

		}
		log.info("Number of reads processed [" + counter + "]");

		if (this.CONTIG_REPORT!=null) writeContigReport(this.CONTIG_REPORT, contigResults);
		if (readWriter!=null) readWriter.close();
		
		// if (this.GENE_REPORT!=null) writeGeneReport (this.GENE_REPORT, geneResults);
		return 0;
	}
	
	
	private void writeReadHeader(List<String> tagNames, PrintStream writer) {
		List<String> header = new ArrayList<>();
		header.add("INPUT_1="+this.INPUT_1.toString());
		header.add("INPUT_2="+this.INPUT_2.toString());
		header.add("READ_QUALITY="+this.READ_QUALITY);
		header.add("TRIM_CONTIG_STRING="+this.TRIM_CONTIG_STRING);
		String h = StringUtils.join(header, "\t");		
		writer.print("#");
		writer.println(h);

		List<String> columns = new ArrayList<String>(Arrays.asList("READ_NAME", "CONTIG_1", "POS_1", "MQ_1", "CONTIG_2", "POS_2", "MQ_2", "PRIMARY_MAPPING"));
		for (String tag: tagNames) {
			columns.add(tag+"_1");
			columns.add(tag+"_2");
		}
		writer.println(StringUtils.join(columns, "\t"));
	}
	
	private void writeReadReportLine(List<SAMRecord> r1List, List<SAMRecord> r2List, List<String> tagNames, PrintStream writer) {
		// only the 2nd read can be a multimapper.
		if (!validateReadSetSize(r1List, r2List)) return;
		
		SAMRecord r1 = r1List.get(0);
		
		String r1Contig=r1.getContig().replaceAll(this.TRIM_CONTIG_STRING, "");
		int r1Pos=r1.getAlignmentStart();
		
		for (SAMRecord r: r2List) {
			List<String> line = new ArrayList<>();
			Collections.addAll(line, r1.getReadName(), r1Contig, Integer.toString(r1Pos), Integer.toString(r1.getMappingQuality()));
			String r2Contig = r.getContig();
			if (r2Contig!=null) r2Contig.replaceAll(this.TRIM_CONTIG_STRING, "");
			if (r2Contig==null) r2Contig="NA";
			int r2Pos=r.getAlignmentStart();
			// if the contig and position are the same, don't emit the read
			if (r1Contig.equals(r2Contig) && r1Pos==r2Pos) 
				continue;
			
			Collections.addAll(line, r2Contig, Integer.toString(r2Pos), Integer.toString(r.getMappingQuality()), Boolean.toString(!r.isSecondaryAlignment()));
			for (String tagName: tagNames) {
				line.add(getTagValueAsString(r1, tagName));
				line.add(getTagValueAsString(r, tagName));
			}
			writer.println(StringUtils.join(line, "\t"));
		}		
	}
	
	private String getTagValueAsString (SAMRecord r, String tagName) {
		Object v = r.getAttribute(tagName);
		if (v==null) return ("NA");
		return v.toString();
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
		
		// handle nulls in contig names.  List of contigs may be empty
		Collection<String> contigsNew = r2.stream().map(x-> x.getContig()).filter(x->x!=null).map(x -> x.replaceAll(this.TRIM_CONTIG_STRING, "")).map(x-> stringInterner.intern(x)).collect(Collectors.toList());
		
		// Collection<String> contigsNew = r2.stream().map(x-> x.getContig()).map(x -> x.replaceAll(this.TRIM_CONTIG_STRING, "")).map(x-> stringInterner.intern(x)).collect(Collectors.toList());
		List<String> contigsNewUnique = new ArrayList<>(new TreeSet<>(contigsNew));
		// if the size is 0 re-encode as NA
		if (contigsNewUnique.size()==0) {
			contigsNewUnique= Collections.singletonList("NA");
		}
		ContigResult r = new ContigResult(contigOne, contigsNewUnique, contigsNewUnique.size()==1);
		return (r);
	}

	private PeekableIterator<List<SAMRecord>> getReadIterator (final List<File> bamFile, final Integer readQuality, boolean removedUnmappedReads) {
        Iterator<SAMRecord> iter = getQueryNameSortedData(bamFile);
		// filter out unmapped reads.
		if (removedUnmappedReads) iter = new UnmappedReadFilter(iter);
		// optionally, filter out reads below a map quality threshold.
		if (readQuality!=null) iter = new MapQualityFilteredIterator(iter, readQuality, false).iterator();
		final GroupingIterator<SAMRecord> groupingIterator = new GroupingIterator<>(iter, READ_NAME_COMPARATOR);
		PeekableIterator<List<SAMRecord>> peekable = new PeekableIterator<>(groupingIterator);
		return peekable;
	}

	private Iterator<SAMRecord> getQueryNameSortedData (final List<File> bamFiles) {
		SamReaderFactory factory= SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE);
		
		SamHeaderAndIterator headerIterator= SamFileMergeUtil.mergeInputs(bamFiles, false, factory);
				
		if (headerIterator.header.getSortOrder().equals(SortOrder.queryname))
			return headerIterator.iterator;
		log.info("Input SAM/BAM not in queryname order, sorting...");
        final ProgressLogger progressLogger = new ProgressLogger(log, 1000000, "Sorting reads in query name order");
        final CloseableIterator<SAMRecord> result = SamRecordSortingIteratorFactory.create(headerIterator.header, headerIterator.iterator, READ_NAME_COMPARATOR, progressLogger);
        log.info("Sorting finished.");
        return result;
	}

	/*
	private void writeGeneReport (final File outFile, final Map<String, GeneResult> geneResults) {
		PrintStream writer = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(outFile));
		List<String> header = new ArrayList<>();
		header.add("INPUT_1="+this.INPUT_1.toString());
		header.add("INPUT_2="+this.INPUT_2.toString());
		header.add("READ_QUALITY="+this.READ_QUALITY);
		// header.add("GENE_EXON_TAG="+this.GENE_EXON_TAG);
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
	*/
	
	private void writeContigReport (final File outFile, final ObjectCounter<ContigResult> contigResults) {
		PrintStream writer = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(outFile));
		List<String> header = new ArrayList<>();
		header.add("INPUT_1="+this.INPUT_1.toString());
		header.add("INPUT_2="+this.INPUT_2.toString());
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
    	private final Log log = Log.getInstance(UnmappedReadFilter.class);
    	
    	public UnmappedReadFilter(final Iterator<SAMRecord> underlyingIterator) {
            super(underlyingIterator);
    	}
        @Override
        public boolean filterOut(final SAMRecord r) {
        	return r.getReadUnmappedFlag();
        }
		@Override
		public void logFilterResults() {
			String msg = String.format("Records pass [%d] records fail [%d] ",this.getRecordsPassed(), this.getRecordsFailed());  
			log.info(msg);											
		}
    }


	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new CompareDropSeqAlignments().instanceMain(args));
	}
}


