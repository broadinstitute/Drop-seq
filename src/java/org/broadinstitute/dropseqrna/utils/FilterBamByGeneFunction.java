package org.broadinstitute.dropseqrna.utils;

import java.io.File;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.GeneFunctionCommandLineBase;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.readiterators.GeneFunctionFilteringIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.MissingTagFilteringIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;
import org.broadinstitute.dropseqrna.utils.readiterators.TagValueFilteringIterator;

import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import picard.annotation.LocusFunction;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        summary = "Filter reads by gene annotations",
        oneLineSummary = "Filter reads by gene annotations",
        programGroup = DropSeq.class
)
public class FilterBamByGeneFunction extends GeneFunctionCommandLineBase {

	private final Log log = Log.getInstance(FilterBamByGeneFunction.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze. This argument can accept wildcards, or a file with the suffix .bam_list that contains the locations of multiple BAM files", minElements = 1)
	public List<File> INPUT;
	
	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output report")
	public File OUTPUT;

	@Argument(doc="Emit reads that have the cell barcodes in this file.  The file has 1 column with no header.", optional=true)
	public File CELL_BC_FILE=null;
	
	@Argument(doc="The cell barcode tag.  If there are no reads with this tag, the program will assume that all reads belong to the same cell and process in single sample mode.")
	public String CELL_BARCODE_TAG="XC";
	
	@Override
	protected int doWork() {
		
		this.INPUT = FileListParsingUtils.expandFileList(INPUT);		
		IOUtil.assertFileIsWritable(OUTPUT);
		
		Collection<String> cellBarcodes = getCellBarcodes();
		
		SamHeaderAndIterator headerAndIter= SamFileMergeUtil.mergeInputs(this.INPUT, false, SamReaderFactory.makeDefault());
		SAMFileWriter out = new SAMFileWriterFactory().makeSAMOrBAMWriter(headerAndIter.header, true, OUTPUT);

		Iterator<SAMRecord> readIter = getIterator(headerAndIter, cellBarcodes, this.GENE_NAME_TAG, this.GENE_STRAND_TAG, 
				this.GENE_FUNCTION_TAG, this.STRAND_STRATEGY, this.LOCUS_FUNCTION_LIST, this.CELL_BARCODE_TAG);
		
		readIter.forEachRemaining(x-> out.addAlignment(x));
		out.close();		
		log.info("Finished");
		
		return 0;
	}
	
	private Iterator<SAMRecord> getIterator (final SamHeaderAndIterator headerAndIter, Collection<String> cellBarcodes, final String geneNameTag,
            final String geneStrandTag, final String geneFunctionTag, final StrandStrategy strandStrategy,
            final Collection <LocusFunction> acceptedLociFunctions, final String cellBarcodeTag) {
		

		final ProgressLogger progressLogger = new ProgressLogger(log);
		Iterator<SAMRecord> filteringIterator = new ProgressLoggingIterator(headerAndIter.iterator, progressLogger);
				
        // Filter records before sorting, to reduce I/O
		filteringIterator = new MissingTagFilteringIterator(filteringIterator, cellBarcodeTag, geneNameTag);

		// Filter reads on if the read contains a cell barcode, if cell barcodes have been specified.
		if (cellBarcodes != null) {
			filteringIterator = new TagValueFilteringIterator<>(filteringIterator, cellBarcodeTag, cellBarcodes);
		}
		// Filter/assign reads based on functional annotations
		filteringIterator = new GeneFunctionFilteringIterator(filteringIterator, geneNameTag, geneStrandTag, geneFunctionTag, strandStrategy, acceptedLociFunctions); 				
		return filteringIterator;
	}
	
	private Collection<String> getCellBarcodes () {
		if (this.CELL_BC_FILE==null) return (null);
		Collection<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(this.CELL_BC_FILE);
		log.info("Found " + cellBarcodes.size()+ " cell barcodes in file");
		return (cellBarcodes);
	}
	
	
	
	public static void main(final String[] args) {
		System.exit(new FilterBamByGeneFunction().instanceMain(args));
	}
}
