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
package org.broadinstitute.dropseqrna.metagene;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.TranscriptomeException;
import org.broadinstitute.dropseqrna.barnyard.GeneFunctionCommandLineBase;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.SamHeaderUtil;
import org.broadinstitute.dropseqrna.utils.SamWriterSink;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import picard.annotation.LocusFunction;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.File;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

@CommandLineProgramProperties(
        summary = "A special case tagger.  Using primary and secondary alignments, try to discover genes that have high homology, where reads map to multiple genes consistently."
        		+ "This takes advantage of UMI barcodes, by looking for UMIs that consistently map to multiple genes.",
        oneLineSummary = "Discovery and tag meta-genes",
        programGroup = DropSeq.class
)

public class DiscoverMetaGenes extends GeneFunctionCommandLineBase {

	private final Log log = Log.getInstance(DiscoverMetaGenes.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze")
	public File INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output BAM, updated with new gene name/strand/function tags for UMIs "
			+ "that are best explained as metagenes.", optional=true)
	public File OUTPUT;

	@Argument (doc="A report on which genes are associated to each other.  Only generated if KNOWN_META_GENE_FILE is not set.", optional=true)
	public File REPORT;

	@Argument(doc="A list of known meta genes.  If this list is provided, new meta genes will not be discovered in the data set, and only meta genes in this list will be tagged.  "
			+ "Formatted as a single column colon seperated list of genes, 1 line per metagene.  This supports partial matching, where the known meta genes are supersets of the "
			+ "discovered metagenes.  For example, if a known metagene is A:B:C and in this data set A:B or A:C is discovered, the results are tagged as A:B:C. Note that if "
			+ "two known metagenes exist and there is an exact and inexact match [A:B and A:B:C in the known list] the exact match is preferred. "
			+ "This can be a previous REPORT file.", optional=true)
	public File KNOWN_META_GENE_FILE;

	@Argument(doc="Should all reads be written to the output BAM, or only the low MQ reads that support metagenes?  When set to true, the output BAM has the same set of reads as the input."
			+ "This process is significantly slower when this flag is set to true.  When set to false, a much smaller BAM file is produced containing ONLY the metagene reads. Note that when "
			+ "this is set to true, uniquely mapped reads on the selected cell barcodes may have modified gene name/strand/functions, as the those tags are filtered by the accepted locus functions "
			+ " and strand strategy.")
	public boolean WRITE_ALL_READS=false;

	@Argument (doc="A list of cell barcodes to exclusively look at.  Use the list of cell barcodes you believe are true cells to best find paired genes during discovery.", optional=true)
	public File CELL_BC_FILE;

	@Argument(doc="The cell barcode tag.  Only used in conjunction with CELL_BC_FILE", optional=true)
	public String CELL_BARCODE_TAG="XC";

	@Argument(doc="The molecular barcode tag.", optional=true)
	public String MOLECULAR_BARCODE_TAG="XM";

	@Argument (doc="The metagene name tag to add to records")
	public String METAGENE_NAME="mn";
	
	@Argument (doc="The metagene strand tag to add to records")
	public String METAGENE_STRAND="ms";
	
	@Argument (doc="The metagene function tag to add to records")
	public String METAGENE_FUNCTION="mf";
		
	@Argument(doc="The minimum map quality of the read to get a metatag.  Using STAR, use a MQ=3 for reads that map to two locations, and a MQ=1 for reads that map to three locations.")
	public Integer MIN_READ_MQ=2;

	@Argument (doc="The minimum mapping quality for a read to be considered uniquely mapped.  Using STAR, this is any number greater than 3.")
	public Integer MAP_QUALITY_UNIQUE=10;

	@Argument(doc="Only consider reads that are within <MAXIMUM_READ_EDIT_DISTANCE> of the place they are mapped to on the reference. If left unset, no reads are filtered by edit distance."
			+ "This filter only applies to model discovery.  All reads regardless of edit distance to the reference will be tagged if they are part of a valid metagene model.", optional=true)
	public Integer MAXIMUM_READ_EDIT_DISTANCE=0;
	
	@Argument(doc="The ratio of unambiguous UMIs to uniquely mapped UMIs for a metagene to be included.  "
			+ "This parameter is only considered if there's no pre-existing set of meta genes [KNOWN_META_GENE_FILE=null].  "
			+ "The default setting filters a large number of suprious alignments and false positives from the data.  Seeting this rate closer to one will further "
			+ "reduce the number of metagenes discovered.")
	public Double META_GENE_RATIO=0.01d;

    @Argument(doc="Should single gene UMI counts be written out as separate lines into the REPORT")
    public boolean WRITE_SINGLE_GENES = false;

    private final String BAM_TAG_DELIMITER=",";
		
	private static final String DELIMITER=":";
	private static final int progressLogInterval=1000;

	@Override
	protected int doWork() {

		IOUtil.assertFileIsReadable(this.INPUT);
		if (this.REPORT!=null) IOUtil.assertFileIsWritable(this.REPORT);
		if (this.OUTPUT!=null) IOUtil.assertFileIsWritable(this.OUTPUT);

		List<String> cellBarcodes=null;
		if (this.CELL_BC_FILE!=null) {
			IOUtil.assertFileIsReadable(this.CELL_BC_FILE);
			cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(this.CELL_BC_FILE);
		}

		// returns null if this.KNOWN_META_GENE_FILE is null.
		Set<MetaGene> discoveredMetaGenes = parseMetaGeneFile(this.KNOWN_META_GENE_FILE);

		// if no metagenes are supplied, discover them.
		if (discoveredMetaGenes==null)
			discoveredMetaGenes= umiLevelMetaGeneDiscovery(cellBarcodes, this.REPORT);

		// tag BAM if output supplied.
		if (this.OUTPUT!=null && discoveredMetaGenes!=null)
			umiLevelMetaGeneTagging(cellBarcodes, discoveredMetaGenes);

		return 0;
	}

	private Set<MetaGene> parseMetaGeneFile (final File input) {
		if (input==null) return null;
		UMIMetaGeneAggregation aggregator = DiscoverMetaGenes.readReport(input);
		Set<MetaGene> discoveredMetaGenes = aggregator.getMetaGeneList(this.META_GENE_RATIO);
		log.info("Discovered " + discoveredMetaGenes.size()+ " meta genes");
		return discoveredMetaGenes;
	}
	
	/**
	 * For each cell barcode + molecular barcode, determine which genes should be meta genes.
	 * @param cellBarcodes
	 * @return
	 */
	public Set<MetaGene> umiLevelMetaGeneDiscovery(final List<String> cellBarcodes, final File report) {
		SamHeaderAndIterator headerAndIterator = SamFileMergeUtil.mergeInputs(Collections.singletonList(this.INPUT), false);
		SamHeaderUtil.addPgRecord(headerAndIterator.header, this);

		UMIMetaGeneAggregation aggregator = new UMIMetaGeneAggregation();
		log.info("Searching for meta genes");
		
		// filters out single reads that map to multiple genes on the same strand
		UMIMetaGeneCollectionIterator iter = new UMIMetaGeneCollectionIterator(headerAndIterator, this.GENE_NAME_TAG, this.GENE_STRAND_TAG, this.GENE_FUNCTION_TAG, 
				this.STRAND_STRATEGY, this.LOCUS_FUNCTION_LIST, this.FUNCTIONAL_STRATEGY, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG,
				this.MIN_READ_MQ, this.MAP_QUALITY_UNIQUE, this.MAXIMUM_READ_EDIT_DISTANCE, cellBarcodes, null, false);
		
		int totalCellBarcodes=0;
		if (cellBarcodes!=null) totalCellBarcodes=cellBarcodes.size();
		ProgressLoggingIterator pl = new ProgressLoggingIterator(iter, log, progressLogInterval, totalCellBarcodes);

		while (pl.hasNext()) {
			UMIMetaGeneCollection c = pl.next();
			// you only add data to the aggregator during discovery if there's informative data to add.
			if (c.getNumInformativeReads()>0)
				aggregator.add(c, false);
		}
		pl.close();		
		if (report!=null) writeReport (aggregator, report, WRITE_SINGLE_GENES);

		Set<MetaGene> discoveredMetaGenes = aggregator.getMetaGeneList(this.META_GENE_RATIO);
		log.info("Discovered " + discoveredMetaGenes.size()+ " meta genes");
		return discoveredMetaGenes;
	}

	public int umiLevelMetaGeneTagging (final List<String> cellBarcodes, final Set<MetaGene> metaGenes) {
		SamHeaderAndIterator headerAndIterator = SamFileMergeUtil.mergeInputs(Collections.singletonList(this.INPUT), false);		
		SamHeaderUtil.addPgRecord(headerAndIterator.header, this);
		// set up the output writer and the sink.
		SAMFileWriter writer=null;
		SamWriterSink sink=null;	
		if (this.OUTPUT!=null)
			writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(headerAndIterator.header, false, OUTPUT);	
		if (this.OUTPUT!=null & this.WRITE_ALL_READS)
			sink= new SamWriterSink(writer);
		 		
		UMIMetaGeneCollectionIterator iter = new UMIMetaGeneCollectionIterator(headerAndIterator, this.GENE_NAME_TAG, this.GENE_STRAND_TAG, this.GENE_FUNCTION_TAG, 
				this.STRAND_STRATEGY, this.LOCUS_FUNCTION_LIST, this.FUNCTIONAL_STRATEGY, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG,
				this.MIN_READ_MQ, this.MAP_QUALITY_UNIQUE, null, cellBarcodes, sink, false);
				
		int totalCellBarcodes=0;
		if (cellBarcodes!=null) totalCellBarcodes=cellBarcodes.size();
		ProgressLoggingIterator pl = new ProgressLoggingIterator(iter, log, progressLogInterval, totalCellBarcodes);

		log.info("Tagging BAM file");

		while (pl.hasNext()) {
			UMIMetaGeneCollection c = pl.next();
			writeMetaGenesToBAM(c, metaGenes, writer, this.WRITE_ALL_READS);
		}		
		pl.close();
		if (writer!=null) {
			writer.close();
		}
		return 0;
	}


	private void writeMetaGenesToBAM (final UMIMetaGeneCollection c, final Collection <MetaGene> approvedMetaGenes, final SAMFileWriter writer , final boolean writeAllReads) {
		// For unambiguous meta genes, tag one of the reads with the new metagene tag.
		// don't need to remove reads, as other reads with the same queryname will not get the metagene tag and will not be processed.
		// should explicitly remove any existing metagene tag from all other reads.
		// when tagging reads, if you're only writing out the meta gene reads, you can pre-filter the data.  So if write all reads is true, then pre-filter equals false.

		// this runs on only the selected cell barcodes.
		Collection <ReadGroupResult> rgr =  c.getResultsPerRead();
		// filter to read groups that are unambiguous and write them out.
		List<ReadGroupResult> metaGeneGroups= rgr.stream().filter(x -> x.getType()==ReadGroupResult.MetaGeneTypeEnum.UNAMBIGIOUS).toList();
		// process metaGenes.  MetaGene has to be on the "approved" list to be tagged.
		for (ReadGroupResult r: metaGeneGroups) {

            // pick the first read and add the metagene tag if it's in the approved list.
            for (SAMRecord rec : r.getReads()) {
                MetaGene m = r.getMetaGene();
                m = getMetaGeneInKnownModel(m, approvedMetaGenes);
                if (m != null)
                    rec = tagMetaGene(rec, m);
                writer.addAlignment(rec);
            }
		}

		// all other read groups that are not meta genes.
		if (writeAllReads) {
			List<ReadGroupResult> other= rgr.stream().filter(x -> x.getType()!=ReadGroupResult.MetaGeneTypeEnum.UNAMBIGIOUS).toList();
			for (ReadGroupResult r: other)
				r.getReads().stream().forEach(writer::addAlignment);
		}

	}
	
	/**
	 * This handles partial matching, as deferred to in MetaGene.partialMatch
	 * 
	 * @param m The metagene result to test
	 * @param approvedMetaGenes A collection of metagenes that are valid
	 * @return The MetaGene that is valid for this result (which may not be the input metagene if the input is a subset of the result), or null if no valid metagene is found.
	 */
	MetaGene getMetaGeneInKnownModel (MetaGene m, Collection<MetaGene> approvedMetaGenes) {
		if (approvedMetaGenes.contains(m))
			return m;
		for (MetaGene a: approvedMetaGenes) {
			if (a.partialMatch(m)) return a;
		}
		return null;
	}
	
	/**
	 * Tag the read with the meta gene tags for the name, strand, and function.  
	 * Since these function as a tuple, need to have the same number of copies of each record.
	 * 
	 * @param rec The SAM record to tag
	 * @param m The meta gene to tag on the SAM record.
	 * @return The tagged SAM record
	 */
	private SAMRecord tagMetaGene (SAMRecord rec, MetaGene m) {
		char sep = DiscoverMetaGenes.DELIMITER.charAt(0);
		
		// get the list of locus functions for the read.  This can be more than 1, like CODING,UTR
		// this means that all 3 tags need to be the same length!
		// split up the function into the list of values.
		List<LocusFunction> funcs = getLocusFunctionsForMetaGene (m, rec);

		// have as many metagene tags as locus function tags.
		String metaGeneName = StringUtils.join(Collections.nCopies(funcs.size(), m.getMetaGeneName(sep)), this.BAM_TAG_DELIMITER);
		rec.setAttribute(this.METAGENE_NAME, metaGeneName);
		
		// Have as many strand tags as locus function tags.					
		String metaGeneStrand = StringUtils.join(Collections.nCopies(funcs.size(), rec.getStringAttribute(this.GENE_STRAND_TAG)), this.BAM_TAG_DELIMITER);
		rec.setAttribute(this.METAGENE_STRAND, metaGeneStrand);
		
		// add the locus function string for the meta gene.
		String locusFunctionString = StringUtils.join(funcs.stream().map(x->x.toString()).collect(Collectors.toList()), this.BAM_TAG_DELIMITER);
		rec.setAttribute(this.METAGENE_FUNCTION, locusFunctionString); 
		return (rec);
	}
	
	/**
	 * When constructing the locus function for the output meta gene, it should be any locus function(s) for the genes the read participates in that overlap the metagene genes.
	 * The read might overlap both the meta gene, as well as some ancillary read (for example on the opposite strand), and those "extra" locus functions should be excluded.
	 * @param m The meta gene this read belongs to.
	 * @param rec The SAM record to process
	 * @return The list of locus functions for this meta gene.
	 */
	private List<LocusFunction> getLocusFunctionsForMetaGene (MetaGene m, SAMRecord rec ) {				
		String function = rec.getStringAttribute(this.GENE_FUNCTION_TAG);
		List<LocusFunction> funcs = Arrays.stream(function.split(BAM_TAG_DELIMITER)).map(LocusFunction::valueOf).toList();
		String []  genes = rec.getStringAttribute(this.GENE_NAME_TAG).split(this.BAM_TAG_DELIMITER);
		
		List<LocusFunction> result = new ArrayList<>();
		
		for (int i=0; i<genes.length; i++) {
			String gene = genes[i];
			LocusFunction func = funcs.get(i);
			if (m.getGeneNames().contains(gene))
				result.add(func);			
		}		
		return result;
	}

    static void writeReportHeader(PrintStream out) {
        String [] line={"GENE_NAME", "GENE_UNIQUE_UMI", "TOTAL_UNIQUE_UMI", "METAGENE_UMI", "AMBIGUOUS_UMI", "FRAC_METAGENE_UMI"};
        String h = StringUtils.join(line, "\t");
        out.println(h);
    }

    public static UMIMetaGeneAggregation readReport(final File reportFile) {
        UMIMetaGeneAggregation aggregator = new UMIMetaGeneAggregation();

        TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(reportFile);
        for (TabbedTextFileWithHeaderParser.Row row : parser) {
            String[] geneNames = StringUtils.split(row.getField("GENE_NAME"), DELIMITER);
            int[] geneUniqueUMIs = Stream.of(StringUtils.split(row.getField("GENE_UNIQUE_UMI"), DELIMITER)).mapToInt(Integer::parseInt).toArray();
            if (geneUniqueUMIs.length != geneNames.length) {
                throw new TranscriptomeException("Error parsing MetaGenes report " +
                    reportFile.getAbsolutePath() + ": " +
                    "geneUniqueUMIs.length [" + geneUniqueUMIs.length + "] is not equal to geneNames.length [" + geneNames.length + "]" +
                    ". GENE_NAME='" + row.getField("GENE_NAME") + "', GENE_UNIQUE_UMI='" + row.getField("GENE_UNIQUE_UMI") + "'");
            }

            ObjectCounter<String> uniqueGenes = aggregator.getUniqueGenes();
            for (int idx=0; idx<geneNames.length; idx++) {
                uniqueGenes.setCount(geneNames[idx], geneUniqueUMIs[idx]);
            }

            // If this is a single gene, don't add it to the meta gene counters
            int metageneUMI = row.getIntegerField("METAGENE_UMI");
            int ambiguousUMI = row.getIntegerField("AMBIGUOUS_UMI");
            if (metageneUMI > 0 || ambiguousUMI > 0) {
                MetaGene metagene = new MetaGene(Arrays.asList(geneNames));
                aggregator.getUnAmbiguousGenes().incrementByCount(metagene, metageneUMI);
                aggregator.getAmbiguousGenes().incrementByCount(metagene, ambiguousUMI);
            }
        }
        parser.close();

        return aggregator;
    }

    static void writeReportLine(String metaGeneName,
                                String uniqueCounts,
                                int totalUniqueCounts,
                                int unambigCount,
                                int ambigCount,
                                double ratio,
                                PrintStream out,
                                DecimalFormat percentageFormat) {
        String [] b ={metaGeneName, uniqueCounts, Integer.toString(totalUniqueCounts), Integer.toString(unambigCount), Integer.toString(ambigCount), percentageFormat.format(ratio)};
        String r = StringUtils.join(b, "\t");
        out.println(r);
    }

	public static void writeReport (final UMIMetaGeneAggregation aggregator, final File outFile, boolean writeSingleGenes) {
		DecimalFormat percentageFormat = new DecimalFormat("###.##");
		PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(outFile));
		char sep = DELIMITER.charAt(0);

        writeReportHeader(out);

        ObjectCounter<String> unique = aggregator.getUniqueGenes();
        ObjectCounter<MetaGene> unambig = aggregator.getUnAmbiguousGenes();
        ObjectCounter<MetaGene> ambig = aggregator.getAmbiguousGenes();
        List<MetaGene> metaGenes = unambig.getKeysOrderedByCount(true);

        if (writeSingleGenes) {
            Set<String> genesSeen = new HashSet<>();
            for (MetaGene metaGene : metaGenes) {
                genesSeen.addAll(metaGene.getGeneNames());
            }
            for (String gene: unique.getKeysOrderedByCount(true)) {
                if (!genesSeen.contains(gene)) {
                    final int umiCount = unique.getCountForKey(gene);
                    writeReportLine(gene, Integer.toString(umiCount), umiCount, 0, 0, 0d, out, percentageFormat);
                }
            }
        }

        for (MetaGene mg: metaGenes) {
        	Set<String> genes = mg.getGeneNames();
        	String metaGeneName = mg.getMetaGeneName(sep);
        	String [] uniqueCountsString = genes.stream().mapToInt(unique::getCountForKey).mapToObj(Integer::toString).toArray(String[]::new);
        	int totalUniqueCounts = genes.stream().mapToInt(unique::getCountForKey).sum();
        	String uniqueCounts= StringUtils.join(uniqueCountsString, sep);
        	int ambigCount = ambig.getCountForKey(mg);
        	int unambigCount = unambig.getCountForKey(mg);
        	double ratio = (double) unambigCount / ((double) totalUniqueCounts + (double) unambigCount);
        	if (Double.isInfinite(ratio)) ratio= unambigCount; // if the number of unique reads is 0, treat it as 1 for this so you don't get infinite ratios.
            writeReportLine(metaGeneName, uniqueCounts, totalUniqueCounts, unambigCount, ambigCount, ratio, out, percentageFormat);
        }
        out.close();
	}

	private class ProgressLoggingIterator extends FilteredIterator<UMIMetaGeneCollection> {
		private final Log log;
		private final Set<String> cells;
		private String currentCell;
		private final int interval;
		private final int numCellsExpected;
		public ProgressLoggingIterator (final Iterator<UMIMetaGeneCollection> underlyingIterator, final Log log, final int interval, final int numCellsExpected) {
			super(underlyingIterator);
			this.log=log;
			this.cells=new HashSet<>();
			this.currentCell="";
			this.interval=interval;
			this.numCellsExpected=numCellsExpected;
		}

		@Override
		public boolean filterOut(final UMIMetaGeneCollection rec) {
			String newCell=rec.getCellBarcode();
			if (newCell==null)
				log.info("STOP");
			if (!newCell.equals(currentCell)) {
				cells.add(newCell);
				if (cells.size()%this.interval==0 || cells.size()==this.numCellsExpected) log.info("Cells processed [" + cells.size() + "/" +  numCellsExpected+"]");
			}
			currentCell=newCell;
			return false;
		}
		
	}

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new DiscoverMetaGenes().instanceMain(args));
	}
}


