package org.broadinstitute.dropseqrna.annotation.functionaldata.disambiguate;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalData;
import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalDataProcessorStrategy;
import org.broadinstitute.dropseqrna.barnyard.*;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.MultiComparator;
import org.broadinstitute.dropseqrna.utils.StringTagComparator;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.*;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Starsolo and DropSeq process some functional annotations differently.  Can these ambiguities be resolved" +
                "by looking at many reads from the same UMIs?  In particular, for reads that are on the opposite strand of a coding" +
                " annotation and on the same strand as an intron, STARsolo interprets those reads as antisense coding, while DropSeq interprets " +
                "those reads as intronic.  Perhaps other reads coming from the same UMI can disambiguate these assignments.",
        oneLineSummary = "Disambiguate UMI functional annotations",
        programGroup = DropSeq.class
)
public class OptimusDropSeqLocusFunctionComparison extends GeneFunctionCommandLineBase {

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze. This argument can accept wildcards, or a file with the suffix .bam_list that contains the locations of multiple BAM files", minElements = 1)
    public List<File> INPUT;

    @Argument(doc="Contains one row for each cell/UMI that has an ambiguous read, and the status of how the UMI is resolved." +
            "Cell/UMI read collections that are unambiguous are not emitted.", optional=true)
    public File OUTPUT=null;

    @Argument(doc="A summary all cell/UMI results.", optional=true)
    public File SUMMARY=null;

    @Argument(doc="File containing a list of cell barcodes to process.  If not provided, process all cell barcodes", optional=true)
    public File CELL_BC_FILE=null;

    @Argument(doc="The cell barcode tag.  If there are no reads with this tag, the program will assume that all reads belong to the same cell and process in single sample mode.")
    public String CELL_BARCODE_TAG="XC";

    @Argument(doc="The molecular barcode tag.")
    public String MOLECULAR_BARCODE_TAG="XM";

    @Argument(doc="The map quality of the read to be included.  The code is not tested (and will no doubt be more complex)" +
            "if non-unique reads are included.")
    public Integer READ_MQ=10;

//    @Argument(doc="Gene Name tag.  Takes on the gene name this read overlaps (if any)")
//    public String GENE_NAME_TAG= DEFAULT_GENE_NAME_TAG;
//
//    @Argument(doc="Gene Strand tag.  For a given gene name <GENE_NAME_TAG>, this is the strand of the gene.")
//    public String GENE_STRAND_TAG= DEFAULT_GENE_STRAND_TAG;
//
//    @Argument(doc="Gene Function tag.  For a given gene name <GENE_NAME_TAG>, this is the function of the gene at this read's position: UTR/CODING/INTRONIC/...")
//    public String GENE_FUNCTION_TAG= DEFAULT_GENE_FUNCTION_TAG;

    private GeneFunctionProcessor gfp;

    private static final Log log = Log.getInstance(OptimusDropSeqLocusFunctionComparison.class);

    @Override
    protected int doWork() {

        ErrorCheckingPrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(this.OUTPUT));
        writeOutputHeader(out);

        // the FunctionalDataProcessorStrategy shouldn't matter here.
        gfp = new GeneFunctionProcessor(this.GENE_NAME_TAG, this.GENE_STRAND_TAG, this.GENE_FUNCTION_TAG,false,
                this.STRAND_STRATEGY, this.LOCUS_FUNCTION_LIST, FunctionalDataProcessorStrategy.STARSOLO);

        Collection <String> cellBarcodes= getCellBarcodes();

        Iterator<List<SAMRecord>> group = getSamRecordIterator(cellBarcodes);

        DisambiguateFunctionalAnnotation dfa = new DisambiguateFunctionalAnnotation();

        // TODO: Gather summary metrics maybe in a metrics object/file?
        while (group.hasNext()) {
            List<SAMRecord> recs= group.next();
            Map<String, List<FunctionalData>> fdMap = getFunctionalAnnotations(recs);
            Map<String, List<FunctionalData>> fdMap2 = getFunctionalAnnotations2(recs);

            String cellBarcode = recs.getFirst().getStringAttribute(this.CELL_BARCODE_TAG);
            String umiBarcode = recs.getFirst().getStringAttribute(this.MOLECULAR_BARCODE_TAG);
            DisambiguationScore score = dfa.run(cellBarcode, umiBarcode, fdMap);
            // do reporting and summarizing
            writeOutputLine(score, out);
        }


        return 0;
    }

    private void writeOutputHeader (PrintStream out) {
        List<String> line = new ArrayList<>();
        line.add("HEADER");
        String h = StringUtils.join(line, "\t");
        out.println(h);
    }

    private void writeOutputLine (DisambiguationScore score, PrintStream out) {
        List<String> line = new ArrayList<>();
        line.add("FOO");
        String h = StringUtils.join(line, "\t");
        out.println(h);
    }

    Map<String, List<FunctionalData>> getFunctionalAnnotations(List<SAMRecord> recs) {
        Map<String, List<FunctionalData>> result = new HashMap<>();
        for (SAMRecord rec: recs) {
            String readName = rec.getReadName();
            List<FunctionalData> fdList = result.get(readName);
            if (fdList==null)
                fdList = new ArrayList<>();
            List<FunctionalData> fd = getFunctionalData(rec);
            fdList.addAll(fd);
            result.put(readName, fdList);
        }
        return (result);
    }


    Map<String, List<FunctionalData>> getFunctionalAnnotations2(List<SAMRecord> recs) {
        return recs.stream()
                .collect(Collectors.groupingBy(SAMRecord::getReadName,
                        Collectors.mapping(this::getFunctionalData, Collectors.flatMapping(List::stream, Collectors.toList()))));
    }

    private List<FunctionalData> getFunctionalData (SAMRecord r) {
        return gfp.getReadFunctions(r, false);
    }


    /**
     * Create an iterator that produces all reads for a cell and molecular barcode.
     * @param cellBarcodes
     * @return
     */
    Iterator<List<SAMRecord>> getSamRecordIterator (Collection <String> cellBarcodes) {
        ProgressLogger prog = new ProgressLogger(log);

        SamHeaderAndIterator headerAndIterator = SamFileMergeUtil.mergeInputs(this.INPUT, false);

        final StringTagComparator cellBarcodeTagComparator = new StringTagComparator(this.CELL_BARCODE_TAG);
        final StringTagComparator umiTagComparator = new StringTagComparator(this.MOLECULAR_BARCODE_TAG);
        final MultiComparator<SAMRecord> multiComparator = new MultiComparator<>(cellBarcodeTagComparator, umiTagComparator);

        // Filter records before sorting, to reduce I/O
        Iterator<SAMRecord> filteringIterator =
                new MissingTagFilteringIterator(headerAndIterator.iterator, this.CELL_BARCODE_TAG, this.GENE_NAME_TAG, this.MOLECULAR_BARCODE_TAG);

        // Filter reads on if the read contains a cell barcode, if cell barcodes have been specified.
        if (cellBarcodes != null) {
            filteringIterator =
                    new TagValueFilteringIterator<>(filteringIterator, this.CELL_BARCODE_TAG, cellBarcodes);
        }

        // filter on read quality.  Let's start with uniquely mapped reads to reduce complexity -
        // multimappers are hard!
        filteringIterator = new MapQualityFilteredIterator(filteringIterator, this.READ_MQ, true);

        CloseableIterator<SAMRecord> sortedAlignmentIterator = SamRecordSortingIteratorFactory.create(
                headerAndIterator.header, filteringIterator, multiComparator, prog);

        Iterator<List<SAMRecord>> group = new GroupingIterator<>(sortedAlignmentIterator, multiComparator);
        return group;
    }


    List<String> getCellBarcodes () {
        List<String> cellBarcodes=new ArrayList<String>();
        if (this.CELL_BC_FILE!=null) {
            cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(this.CELL_BC_FILE);
            log.info("Found " + cellBarcodes.size()+ " cell barcodes in file");
        }
        return (cellBarcodes);
    }
    @Override
    protected String[] customCommandLineValidation() {
        final ArrayList<String> list = new ArrayList<>(1);

        this.INPUT = FileListParsingUtils.expandFileList(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsWritable(SUMMARY);
        if (CELL_BC_FILE!=null)
            IOUtil.assertFileIsReadable(CELL_BC_FILE);
        return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);
    }

        /** Stock main method. */
    public static void main(final String[] args) {
        System.exit(new OptimusDropSeqLocusFunctionComparison().instanceMain(args));
    }

}
