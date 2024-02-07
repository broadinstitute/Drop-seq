package org.broadinstitute.dropseqrna.annotation.functionaldata.disambiguate;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
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
import org.broadinstitute.dropseqrna.annotation.functionaldata.StarSoloFunctionalDataProcessor;
import org.broadinstitute.dropseqrna.barnyard.*;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.*;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.*;
import picard.annotation.LocusFunction;
import picard.cmdline.CommandLineProgram;
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
public class OptimusDropSeqLocusFunctionComparison extends CommandLineProgram {

    @Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze. This argument can accept wildcards, or a file with the suffix .bam_list that contains the locations of multiple BAM files", minElements = 1)
    public List<File> INPUT;

    @Argument(doc="Contains one row for each cell/UMI that has an ambiguous read, and the status of how the UMI is resolved." +
            "Cell/UMI read collections that are unambiguous or unresolvable are not emitted.", optional=false)
    public File OUTPUT=null;

    @Argument(doc="A BAM file that captures reads fro cell/UMI groupings that had discordant tagging between STARsolo " +
            "and Dropseq.  Purely for debugging purposes.", optional = true)
    public File OUT_BAM=null;

    @Argument(doc="A summary all cell/UMI results.", optional=true)
    public File SUMMARY=null;

    @Argument(doc= "For each read, gather the functional data type (coding, intronic, etc) for each read classified by STAR and DropSeq.  Emit a confusion matrix" +
            "of the results", optional = true)
    public File OUT_CONFUSION_MATRIX =null;

    @Argument(doc="File containing a list of cell barcodes to process.  If not provided, process all cell barcodes", optional=true)
    public File CELL_BC_FILE=null;

    @Argument(doc="The cell barcode tag.  If there are no reads with this tag, the program will assume that all reads belong to the same cell and process in single sample mode.")
    public String CELL_BARCODE_TAG="CB";

    @Argument(doc="The molecular barcode tag.")
    public String MOLECULAR_BARCODE_TAG="UB";

    @Argument(doc="The map quality of the read to be included.  The code is not tested (and will no doubt be more complex)" +
            "if non-unique reads are included.")
    public Integer READ_MQ=10;

    private final String VALIDATION_TAG_STAR="VS";
    private final String VALIDATION_TAG_DROPSEQ="VD";

    private String STARSOLO_FUNCTION_TAG="sF";
    private String STARSOLO_GENE_NAME="GN";

    private final String RECORD_SEP="\t";

//    @Argument(doc="Gene Name tag.  Takes on the gene name this read overlaps (if any)")
//    public String GENE_NAME_TAG= DEFAULT_GENE_NAME_TAG;
//
//    @Argument(doc="Gene Strand tag.  For a given gene name <GENE_NAME_TAG>, this is the strand of the gene.")
//    public String GENE_STRAND_TAG= DEFAULT_GENE_STRAND_TAG;
//
//    @Argument(doc="Gene Function tag.  For a given gene name <GENE_NAME_TAG>, this is the function of the gene at this read's position: UTR/CODING/INTRONIC/...")
//    public String GENE_FUNCTION_TAG= DEFAULT_GENE_FUNCTION_TAG;

    private final List<LocusFunction> LOCUS_FUNCTION_LIST = Collections.unmodifiableList(new ArrayList<>(Arrays.asList(LocusFunction.CODING, LocusFunction.UTR, LocusFunction.INTRONIC)));
    private final StrandStrategy STRAND_STRATEGY=GeneFunctionCommandLineBase.DEFAULT_STRAND_STRATEGY;

    private GeneFunctionProcessor gfp;

    private static final Log log = Log.getInstance(OptimusDropSeqLocusFunctionComparison.class);

    private ObjectCounter<Boolean> readTagValidation;

    @Override
    protected int doWork() {
        readTagValidation = new ObjectCounter<>();
        ProgressLogger pl = new ProgressLogger(log);

        // Set up output
        ErrorCheckingPrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(this.OUTPUT));
        writeOutputHeader(out);

        SamHeaderAndIterator headerAndIterator = SamFileMergeUtil.mergeInputs(this.INPUT, false);
        // the writer can be null.
        SAMFileWriter writer = getBamWriter (headerAndIterator.header, this.OUT_BAM);

        // the FunctionalDataProcessorStrategy shouldn't matter here.
        gfp = new GeneFunctionProcessor(GeneFunctionCommandLineBase.DEFAULT_GENE_NAME_TAG, GeneFunctionCommandLineBase.DEFAULT_GENE_STRAND_TAG,
                GeneFunctionCommandLineBase.DEFAULT_GENE_FUNCTION_TAG,false,
                STRAND_STRATEGY, LOCUS_FUNCTION_LIST, FunctionalDataProcessorStrategy.STARSOLO);

        ValidateAnnotations validator = new ValidateAnnotations(gfp);

        Collection <String> cellBarcodes= getCellBarcodes();
        Iterator<List<SAMRecord>> group = getSamRecordIterator(headerAndIterator, cellBarcodes);
        DisambiguateFunctionalAnnotation dfa = new DisambiguateFunctionalAnnotation();
        // Stores the summary results per cell/UMI.
        ObjectCounter<FunctionCategory> resolvedSummaryCounter = new ObjectCounter<>();

        log.info("Scanning cell / UMIs for ambiguous antisense coding / sense intronic scenarios.");
        ProgressLogger plIter = new ProgressLogger(log, 1000000);

        // Gather summary metrics maybe in a metrics object/file?
        ConfusionMatrix<FunctionalData.Type> confusionMatrix = new ConfusionMatrix<>(FunctionalData.Type.class);

        // Operating on all reads for a single cell / UMI.
        outerloop: while (group.hasNext()) {
            List<SAMRecord> recs= group.next();
            recs.forEach(plIter::record);
            // short circuit single read UMIs.
            if (recs.size()<2)
                continue;

            String cellBarcode = recs.getFirst().getStringAttribute(this.CELL_BARCODE_TAG);
            String umiBarcode = recs.getFirst().getStringAttribute(this.MOLECULAR_BARCODE_TAG);

            Map<String, List<FunctionalData>> dropSeqFD = getFunctionalAnnotations(recs);

            // validate that the tags are approximately the same.
            // if even a single read for the cell / UMI does not validate, do not evaluate this result.
            boolean readsValid = validateFunctionalAnnotations (validator, recs, dropSeqFD, confusionMatrix, writer);
            if (!readsValid)
                continue outerloop;

            DisambiguationScore score = dfa.run(cellBarcode, umiBarcode, dropSeqFD);
            // capture the functional category so we can summarize how many of each type we see.
            FunctionCategory cat = score.classify();
            resolvedSummaryCounter.increment(score.classify());
            if (cat==FunctionCategory.RESOLVED_SENSE_INTRONIC || cat==FunctionCategory.RESOLVED_ANTISENSE_CODING) {
                Set<String> contigs = getUniqueContigs(recs);
                ObjectCounter<String> readStrandCounts = getStrandCount (recs);
                writeOutputLine(score, contigs, readStrandCounts, out);
            }
        }

        // capture the differences in classification of read functional data type.
        writeConfusionMatrix(confusionMatrix);

        // write the summary results
        writeSummary (resolvedSummaryCounter);

        // wrap up writing to the output BAM.
        if (writer!=null)
            writer.close();

        // close the per cell/UMI writer.
        out.close();

        return 0;
    }

    private ObjectCounter<String> getStrandCount (List<SAMRecord> recs) {
        ObjectCounter<String> strand = new ObjectCounter<>();
        recs.forEach(x-> strand.increment(Utils.negativeStrandToString(x.getReadNegativeStrandFlag())));
        return strand;
    }



    /**
     * Validate the STARsolo and dropseq annotations to make sure they are the same before further processing.
     * @param validator
     * @param recs
     * @param dropSeqFD
     * @param confusionMatrix
     * @param writer
     * @return
     */
    private boolean validateFunctionalAnnotations (ValidateAnnotations validator, List<SAMRecord> recs, Map<String, List<FunctionalData>> dropSeqFD,
                                                ConfusionMatrix<FunctionalData.Type> confusionMatrix, SAMFileWriter writer) {
        Map<String, FunctionalData> starsoloFD= getStarSoloAnnotations (recs);
        Map<String, ValidationStatus> validationResult = validator.validate (starsoloFD, dropSeqFD, false);
        boolean readsValue=true;
        for (SAMRecord r: recs) {
            ValidationStatus vs = validationResult.get(r.getReadName());
            confusionMatrix.update(vs.getStar().getCategory(), vs.getDropseq().getCategory());
            boolean isValid=vs.isValid();
            readTagValidation.increment(isValid);
            if (!vs.isValid()) {
                r.setAttribute(this.VALIDATION_TAG_STAR, vs.getStar().getCategory().toString());
                r.setAttribute(this.VALIDATION_TAG_DROPSEQ, vs.getDropseq().getCategory().toString());
                writer.addAlignment(r);
                readsValue=false;
            }
        }
        return readsValue;
    }

    private void writeConfusionMatrix (ConfusionMatrix<FunctionalData.Type> m) {
        if (this.OUT_CONFUSION_MATRIX ==null)
            return;

        ErrorCheckingPrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(this.OUT_CONFUSION_MATRIX));
        m.writeFile(out, '\t', "STARsolo in rows, DropSeq in columns.");
        out.close();
    }

    private void writeSummary (ObjectCounter<FunctionCategory> resolvedSummaryCounter) {
        if (this.SUMMARY==null)
            return;

        // prepare the summary file if not null
        ErrorCheckingPrintStream out= new ErrorCheckingPrintStream(IOUtil.openFileForWriting(this.SUMMARY));
        // write the header
        List<String> header = List.of ("CATEGORY", "COUNT");
        String h = StringUtils.join(header, RECORD_SEP);
        out.println(h);

        // write the categories and counts, sorted from largest to smallest.
        List<FunctionCategory> keys = resolvedSummaryCounter.getKeysOrderedByCount(true);
        for (FunctionCategory key: keys) {
            int count = resolvedSummaryCounter.getCountForKey(key);
            List<String> line = List.of (key.toString(), Integer.toString(count));
            String body = StringUtils.join(line, RECORD_SEP);
            out.println(body);
        }
        out.close();;

    }

    /**
     * Starsolo encodes one and only one functional data annotation per read.
     * The tags are written to the SAM file after interpretation.
     * @param recs
     * @return
     */
    private Map<String, FunctionalData> getStarSoloAnnotations (List<SAMRecord> recs) {
        Map<String, FunctionalData> result = new HashMap<>();
        for (SAMRecord r: recs) {
            Object func = r.getAttribute(this.STARSOLO_FUNCTION_TAG);
            int [] funcArray = (int []) func;
            String geneName = r.getStringAttribute(this.STARSOLO_GENE_NAME);
            boolean readNegativeStrand = r.getReadNegativeStrandFlag();
            FunctionalData fd = StarSoloFunctionalDataProcessor.getFunctionalData(geneName, funcArray[0], funcArray[1], readNegativeStrand);
            result.put(r.getReadName(), fd);
        }
        return result;
    }


    private void writeOutputHeader (PrintStream out) {
        List<String> line = Arrays.asList("CELL_BARCODE", "MOLECULAR_BARCODE", "CONTIG", "CLASS", "AMBIGUOUS_ANTISENSE_CODING", "AMBIGUOUS_SENSE_INTRONIC", "ANTISENSE_CODING", "SENSE_INTRONIC", "CODING", "NUM_READS_POSITIVE_STRAND", "NUM_READS_NEGATIVE_STRAND");
        String h = StringUtils.join(line, RECORD_SEP);
        out.println(h);
    }

    private void writeOutputLine (DisambiguationScore score, Set<String> contigs, ObjectCounter<String> readStrandCounts, PrintStream out) {
        String contigString = contigs.toString();
        if (contigs.size()==1)
            contigString=contigs.iterator().next().toString();
        List<String> line = new ArrayList<>();
        line.add(score.getCell());
        line.add(score.getMolecularBarcode());
        line.add(contigString);
        line.add(score.classify().toString());
        line.add(score.getAmbiguousSenseIntronicCount().toString());
        line.add(score.getAmbiguousSenseIntronicCount().toString());
        line.add(score.getAntisenseCodingCount().toString());
        line.add(score.getSenseIntronicCount().toString());
        line.add(score.getSenseCodingCount().toString());

        // line.add(Integer.toString(score.getTotalCount()));
        line.add(Integer.toString(readStrandCounts.getCountForKey("+")));
        line.add(Integer.toString(readStrandCounts.getCountForKey("-")));
        String h = StringUtils.join(line, RECORD_SEP);
        out.println(h);
    }

    /**
     * Get functional data annotations using the DropSeq tags.
     * @param recs The SAMRecords to transform
     * @return A mapping from the gene name to a list of functional annotations.
     */
    Map<String, List<FunctionalData>> getFunctionalAnnotations(List<SAMRecord> recs) {
        return recs.stream()
                .collect(Collectors.groupingBy(SAMRecord::getReadName,
                        Collectors.mapping(this::getFunctionalData, Collectors.flatMapping(List::stream, Collectors.toList()))));
    }

    /**
     * Convert a single SAMRecord with DropSeq functional annotation tags to a list of functional annotation records.
     * @param r The read to convert
     * @return A list of functional annotations.
     */
    private List<FunctionalData> getFunctionalData(SAMRecord r) {
        return gfp.getReadFunctions(r, false);
    }


    /**
     * Create an iterator that produces all reads for a cell and molecular barcode.
     * @param cellBarcodes
     * @return
     */
    Iterator<List<SAMRecord>> getSamRecordIterator (SamHeaderAndIterator headerAndIterator, Collection <String> cellBarcodes) {
        final ProgressLogger logger = new ProgressLogger(log);

        final StringTagComparator cellBarcodeTagComparator = new StringTagComparator(this.CELL_BARCODE_TAG);
        final StringTagComparator umiTagComparator = new StringTagComparator(this.MOLECULAR_BARCODE_TAG);
        final MultiComparator<SAMRecord> multiComparator = new MultiComparator<>(cellBarcodeTagComparator, umiTagComparator);

        // Filter records before sorting, to reduce I/O
        Iterator<SAMRecord> filteringIterator =
                new MissingTagFilteringIterator(headerAndIterator.iterator, this.CELL_BARCODE_TAG,
                        GeneFunctionCommandLineBase.DEFAULT_GENE_NAME_TAG, this.MOLECULAR_BARCODE_TAG);

        // Filter reads on if the read contains a cell barcode, if cell barcodes have been specified.
        filteringIterator = new CellBarcodeFilteringIterator(filteringIterator, this.CELL_BARCODE_TAG, cellBarcodes);

        // filter on read quality.  Let's start with uniquely mapped reads to reduce complexity -multimappers are hard!
        filteringIterator = new MapQualityFilteredIterator(filteringIterator, this.READ_MQ, true);

        CloseableIterator<SAMRecord> sortedAlignmentIterator = SamRecordSortingIteratorFactory.create(
                headerAndIterator.header, filteringIterator, multiComparator, logger);

        Iterator<List<SAMRecord>> group = new GroupingIterator<>(sortedAlignmentIterator, multiComparator);
        return group;
    }


    List<String> getCellBarcodes () {
        if (CELL_BC_FILE!=null) {
            List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(CELL_BC_FILE);
            log.info("Found " + cellBarcodes.size()+ " cell barcodes in file");
            return (cellBarcodes);
        } else {
            log.info("No cell barcodes file provided, using all cell barcodes");
            return Collections.EMPTY_LIST;
        }

    }

    private SAMFileWriter getBamWriter (SAMFileHeader header, File outBAM) {
        if (outBAM==null)
            return null;

        IOUtil.assertFileIsWritable(outBAM);
        SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, outBAM);
        return writer;
    }

    public static Set<String> getUniqueContigs(List<SAMRecord> samRecords) {
        return samRecords.stream()
                .map(SAMRecord::getContig)
                .collect(Collectors.toSet());
    }

    @Override
    protected String[] customCommandLineValidation() {
        final ArrayList<String> list = new ArrayList<>(1);

        this.INPUT = FileListParsingUtils.expandFileList(INPUT);
        if (this.OUTPUT!=null) IOUtil.assertFileIsWritable(OUTPUT);
        if (this.SUMMARY!=null) IOUtil.assertFileIsWritable(SUMMARY);
        if (this.OUT_CONFUSION_MATRIX !=null) IOUtil.assertFileIsWritable(OUT_CONFUSION_MATRIX);
        if (CELL_BC_FILE!=null)
            IOUtil.assertFileIsReadable(CELL_BC_FILE);
        return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);
    }

        /** Stock main method. */
    public static void main(final String[] args) {
        System.exit(new OptimusDropSeqLocusFunctionComparison().instanceMain(args));
    }

}
