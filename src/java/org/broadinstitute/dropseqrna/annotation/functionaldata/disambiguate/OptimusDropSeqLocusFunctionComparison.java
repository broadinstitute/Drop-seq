package org.broadinstitute.dropseqrna.annotation.functionaldata.disambiguate;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.*;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.MultiComparator;
import org.broadinstitute.dropseqrna.utils.StringTagComparator;
import org.broadinstitute.dropseqrna.utils.readiterators.*;
import picard.cmdline.StandardOptionDefinitions;

import java.io.File;
import java.util.*;

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

    @Argument(doc="A summary of each gene that can be disambiguated.", optional=true)
    public File OUTPUT=null;

    @Argument(doc="File containing a list of cell barcodes to process.  If not provided, process all cell barcodes", optional=true)
    public File CELL_BC_FILE=null;

    @Argument(doc="The cell barcode tag.  If there are no reads with this tag, the program will assume that all reads belong to the same cell and process in single sample mode.")
    public String CELL_BARCODE_TAG="XC";

    @Argument(doc="The molecular barcode tag.")
    public String MOLECULAR_BARCODE_TAG="XM";

    @Argument(doc="Gene Name tag.  Takes on the gene name this read overlaps (if any)")
    public String GENE_NAME_TAG= DEFAULT_GENE_NAME_TAG;

    @Argument(doc="Gene Strand tag.  For a given gene name <GENE_NAME_TAG>, this is the strand of the gene.")
    public String GENE_STRAND_TAG= DEFAULT_GENE_STRAND_TAG;

    @Argument(doc="Gene Function tag.  For a given gene name <GENE_NAME_TAG>, this is the function of the gene at this read's position: UTR/CODING/INTRONIC/...")
    public String GENE_FUNCTION_TAG= DEFAULT_GENE_FUNCTION_TAG;

    private static final Log log = Log.getInstance(OptimusDropSeqLocusFunctionComparison.class);

    @Override
    protected int doWork() {
        Collection <String> cellBarcodes= getCellBarcodes();



        return 0;
    }

    GroupingIterator<SAMRecord> getSamRecordIterator (Collection <String> cellBarcodes) {
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

        // filter on read quality!  Let's start with uniquely mapped reads to reduce complexity -
        // multimappers are hard.

        CloseableIterator<SAMRecord> sortedAlignmentIterator = SamRecordSortingIteratorFactory.create(
                headerAndIterator.header, filteringIterator, multiComparator, prog);

        GroupingIterator<SAMRecord> group = new GroupingIterator<>(sortedAlignmentIterator, multiComparator);

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


        return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);

    }

        /** Stock main method. */
    public static void main(final String[] args) {
        System.exit(new OptimusDropSeqLocusFunctionComparison().instanceMain(args));
    }

}
