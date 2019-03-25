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
package org.broadinstitute.dropseqrna.barnyard;

import java.io.File;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.DGELongFormatRecord.CellBarcodeOrderComparator;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeader;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeaderCodec;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeaderLibrary;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.UMIIterator;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.samtools.util.StringUtil;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        summary = "Measures the digital expression of a library.  " +
                "Method: 1) For each gene, find the molecular barcodes on the exons of that gene.  " +
                "2) Determine how many HQ mapped reads are assigned to each barcode.  " +
                "3) Collapse barcodes by edit distance.  " +
                "4) Throw away barcodes with less than threshold # of reads. " +
                "5) Count the number of remaining unique molecular barcodes for the gene." +
                "This program requires a tag for what gene a read is on, a molecular barcode tag, and a exon tag.  The exon and gene tags may not be present on every read." +
                "When filtering the data for a set of barcodes to use, the data is filtered by ONE of the following methods (and if multiple params are filled in, the top one takes precidence):" +
                "1) CELL_BC_FILE, to filter by the some fixed list of cell barcodes" +
                "2) MIN_NUM_GENES_PER_CELL " +
                "3) MIN_NUM_TRANSCRIPTS_PER_CELL " +
                "4) NUM_CORE_BARCODES " +
                "5) MIN_NUM_READS_PER_CELL",
        oneLineSummary = "Calculate Digital Expression",
        programGroup = DropSeq.class
)

public class DigitalExpression extends DGECommandLineBase {

    private static final Log log = Log.getInstance(DigitalExpression.class);

    @Argument(doc="A summary of the digital expression output, containing 3 columns - the cell barcode, the #genes, and the #transcripts.", optional=true)
    public File SUMMARY=null;

    @Argument(doc="Output number of reads instead of number of unique molecular barcodes.", optional=true)
    public boolean OUTPUT_READS_INSTEAD=false;

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="Output file of DGE Matrix.  Genes are in rows, cells in columns.  The first column contains the gene name. This supports zipped formats like gz and bz2.")
    public File OUTPUT;

    @Argument (doc="An alternate output of expression where each row represents a cell, gene, and count of UMIs.  Cell/Gene pairings with 0 UMIs are not emitted.", optional=true)
    public File OUTPUT_LONG_FORMAT;

    @Argument (doc="Output only genes with at least this total expression level, after summing across all cells", optional=true)
    public Integer MIN_SUM_EXPRESSION=null;

    @Argument(shortName = "H", doc="If true, write a header in the DGE file.  If not specified, and UEI is specified, it is set to true.  " +
            "REFERENCE_SEQUENCE only used to write to header.  If it is not present, it is extracted from INPUT header if possible.", optional = true)
    public Boolean OUTPUT_HEADER;

    @Argument(shortName = "UEI", doc="If OUTPUT_HEADER=true, this is required", optional = true)
    public String UNIQUE_EXPERIMENT_ID;

    private boolean OUTPUT_EXPRESSED_GENES_ONLY=false;

    @Override
    /**
     * This is a revision of the original DGE code to implement a more complicated state machine in the main loop and in exchange get rid of the batch system.
     * This change lets you store the counts of molecular barcodes for gene/cell instead of the reads for one gene/cell, which can save huge amounts of memory when
     * cells have millions of reads on housekeeping genes.
     */
    protected int doWork() {

        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        if (OUTPUT_LONG_FORMAT!=null) IOUtil.assertFileIsWritable(OUTPUT_LONG_FORMAT);
        if (OUTPUT_HEADER == null)
			OUTPUT_HEADER = (UNIQUE_EXPERIMENT_ID != null);
        if (this.SUMMARY!=null) IOUtil.assertFileIsWritable(this.SUMMARY);

        if (REFERENCE_SEQUENCE == null && OUTPUT_HEADER) {
            final SAMFileHeader header = SamReaderFactory.makeDefault().open(INPUT).getFileHeader();
            final SAMSequenceRecord sequence = header.getSequence(0);
            if (sequence != null) {
                String uri = sequence.getAttribute(SAMSequenceRecord.URI_TAG);
                if (uri != null) {
                    final String filePrefix = "file:";
                    if (uri.startsWith(filePrefix))
						uri = uri.substring(filePrefix.length());
                    REFERENCE_SEQUENCE = new File(uri);
                }
            }
        }

        // boolean check = new Utils().validateGetCellBarcodeListParams(this.INPUT, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG,
        //		this.GENE_TAG, this.EXON_TAG, this.CELL_BC_FILE, this.READ_MQ, this.MIN_NUM_TRANSCRIPTS_PER_CELL,
        //		this.MIN_NUM_GENES_PER_CELL, this.MIN_NUM_READS_PER_CELL, this.NUM_CORE_BARCODES);


        List<String> cellBarcodes=new BarcodeListRetrieval().getCellBarcodes(this.INPUT, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG,
                this.GENE_NAME_TAG, this.GENE_STRAND_TAG, this.GENE_FUNCTION_TAG, this.STRAND_STRATEGY, this.LOCUS_FUNCTION_LIST,
                this.CELL_BC_FILE, this.READ_MQ, this.MIN_NUM_TRANSCRIPTS_PER_CELL,
                this.MIN_NUM_GENES_PER_CELL, this.MIN_NUM_READS_PER_CELL, this.NUM_CORE_BARCODES, this.EDIT_DISTANCE, this.MIN_BC_READ_THRESHOLD);

        if (cellBarcodes.isEmpty()) {
            log.error("Running digital expression without somehow selecting a set of barcodes to process no longer supported.");
            return (1);
        } else {
            log.info("Calculating digital expression for [" + cellBarcodes.size()+ "] cells.");
            digitalExpression(cellBarcodes);
        }

        return 0;
    }


    @Override
    protected String[] customCommandLineValidation() {
        final String[] superErrors = super.customCommandLineValidation();
        if (OUTPUT_HEADER == null)
			OUTPUT_HEADER = (UNIQUE_EXPERIMENT_ID != null);
        if (UNIQUE_EXPERIMENT_ID != null || !OUTPUT_HEADER)
			return superErrors;
		else {
            final ArrayList<String> list = new ArrayList<>(1);
            if (superErrors != null)
				for (final String msg: superErrors)
					list.add(msg);
            list.add("If OUTPUT_HEADER==true, UNIQUE_EXPERIMENT_ID must be set");
            return list.toArray(new String[list.size()]);
        }
    }

    private void digitalExpression(final List<String> cellBarcodes) {
        PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(OUTPUT));

        if (OUTPUT_HEADER)
			writeDgeHeader(out);

        writeHeader(out, cellBarcodes);
        //TODO should the ambiguous reads handling be a parameter?  It's set to false by default for DGE to get rid of ambiguous gene assignments on reads
        UMIIterator umiIterator = new UMIIterator(SamFileMergeUtil.mergeInputs(Collections.singletonList(this.INPUT), false),
        		GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG, this.STRAND_STRATEGY, this.LOCUS_FUNCTION_LIST,
        		this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, this.READ_MQ, false, cellBarcodes);

        String gene = null;
        Map<String, Integer> transcriptCountMap = new HashMap<>();
        Map<String, Integer> readCountMap = new HashMap<>();
        Map<String, DESummary> summaryMap = initializeSummary(cellBarcodes);

    	// key: cell barcode.  value: counts of UMIs per gene.
        SortingCollection<DGELongFormatRecord> longFormatRecordCollection = null;

        if (this.OUTPUT_LONG_FORMAT!=null)
        	longFormatRecordCollection=makeSortingCollection(cellBarcodes);

        UMICollection batch;
        while ((batch=umiIterator.next())!=null) {
            if (batch==null || batch.isEmpty())
				continue;
            String currentGene = batch.getGeneName();
            // if just starting the loop
            if (gene==null) gene=currentGene;
            // if the gene is the same, you're still gathering expression on that gene.
            if (gene.equals(currentGene)) {

                if (this.RARE_UMI_FILTER_THRESHOLD>0) batch.filterByUMIFrequency(this.RARE_UMI_FILTER_THRESHOLD);
                int molBCCount = batch.getDigitalExpression(this.MIN_BC_READ_THRESHOLD, this.EDIT_DISTANCE, this.OUTPUT_READS_INSTEAD);
                transcriptCountMap.put(batch.getCellBarcode(), molBCCount);
                int readCount = batch.getDigitalExpression(this.MIN_BC_READ_THRESHOLD, this.EDIT_DISTANCE, true);
                readCountMap.put(batch.getCellBarcode(), readCount);

                // if you're gather the long file format, do it here.
                if (longFormatRecordCollection!=null)
                	addLongFormatRecord(longFormatRecordCollection, batch.getCellBarcode(), batch.getGeneName(), molBCCount);
            }
            // you've gathered all the data for the gene, write it out and start on the next.
            if (!gene.equals(currentGene)) {
                writeStats (gene, transcriptCountMap, cellBarcodes, out);
                addToSummary(readCountMap, transcriptCountMap, summaryMap);
                transcriptCountMap.clear();
                // start the next gene
                if (this.RARE_UMI_FILTER_THRESHOLD>0) batch.filterByUMIFrequency(this.RARE_UMI_FILTER_THRESHOLD);
                int molBCCount = batch.getDigitalExpression(this.MIN_BC_READ_THRESHOLD, this.EDIT_DISTANCE, this.OUTPUT_READS_INSTEAD);
                transcriptCountMap.put(batch.getCellBarcode(), molBCCount);
                int readCount = batch.getDigitalExpression(this.MIN_BC_READ_THRESHOLD, this.EDIT_DISTANCE, true);
                readCountMap.put(batch.getCellBarcode(), readCount);
                gene=currentGene;
                // if you're gather the long file format, do it here.
                if (longFormatRecordCollection!=null)
                	addLongFormatRecord(longFormatRecordCollection, batch.getCellBarcode(), batch.getGeneName(), molBCCount);

            }
        }
        // write out remainder
        if (transcriptCountMap.isEmpty()==false) {
            writeStats (gene, transcriptCountMap, cellBarcodes, out);
            addToSummary(readCountMap, transcriptCountMap, summaryMap);
        }
        out.close();
        if (this.SUMMARY!=null)
			writeSummary(summaryMap.values(), this.SUMMARY);

        if (this.OUTPUT_LONG_FORMAT!=null)
			writeLongOutputFormat(longFormatRecordCollection, this.OUTPUT_LONG_FORMAT);
        CloserUtil.close(umiIterator);

    }

    private void addLongFormatRecord (final SortingCollection<DGELongFormatRecord> longFormatRecords, final String cellBarcode, final String gene, final int umiCount) {
    	DGELongFormatRecord r = new DGELongFormatRecord(cellBarcode, gene, umiCount);
    	longFormatRecords.add(r);
    }

    /**
     * Writes the "long" form of the DGE, with each row representing a cell and gene with the number of UMIs.
     * This is only emitted for cell/gene pairs that are non-null.
     * @param outFile The output file to write to.
     */
    private void writeLongOutputFormat(final SortingCollection<DGELongFormatRecord> longFormatRecordCollection, final File outFile) {
    	final CloseableIterator<DGELongFormatRecord> iterator = longFormatRecordCollection.iterator();

    	PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(outFile));
    	String [] header = {"CELL", "GENE", "UMI_COUNT"};
    	out.println(StringUtil.join("\t", header));
    	while (iterator.hasNext()) {
    		DGELongFormatRecord rec = iterator.next();
    		// write out the cell.
    		String [] line = {rec.getCell(), rec.getGene(), Integer.toString(rec.getCount())};
    		out.println(StringUtil.join("\t", line));

    	}
    	out.flush();
    	CloserUtil.close(out);
    }

    private void writeDgeHeader(final PrintStream out) {
        DgeHeader header = new DgeHeader();
        header.setExpressionFormat(DgeHeader.ExpressionFormat.raw);
        DgeHeaderLibrary lib = new DgeHeaderLibrary(UNIQUE_EXPERIMENT_ID);
        if (REFERENCE_SEQUENCE != null)
			lib.setReference(REFERENCE_SEQUENCE.getAbsoluteFile());
        lib.setInput(INPUT.getAbsoluteFile());
        setDgeHeaderLibraryField(lib, "OUTPUT_READS_INSTEAD", OUTPUT_READS_INSTEAD);
        setDgeHeaderLibraryField(lib, "MIN_SUM_EXPRESSION", MIN_SUM_EXPRESSION);
        setDgeHeaderLibraryField(lib, "EDIT_DISTANCE", EDIT_DISTANCE);
        setDgeHeaderLibraryField(lib, "READ_MQ", READ_MQ);
        setDgeHeaderLibraryField(lib, "MIN_BC_READ_THRESHOLD", MIN_BC_READ_THRESHOLD);
        setDgeHeaderLibraryField(lib, "MIN_NUM_READS_PER_CELL", MIN_NUM_READS_PER_CELL);
        setDgeHeaderLibraryField(lib, "MIN_NUM_GENES_PER_CELL", MIN_NUM_GENES_PER_CELL);
        setDgeHeaderLibraryField(lib, "MIN_NUM_TRANSCRIPTS_PER_CELL", MIN_NUM_TRANSCRIPTS_PER_CELL);
        setDgeHeaderLibraryField(lib, "NUM_CORE_BARCODES", NUM_CORE_BARCODES);
        // if the file is null, don't try to get the absolute file, because that's...uncomfortable.
        if (this.CELL_BC_FILE!=null)
			setDgeHeaderLibraryField(lib, "CELL_BC_FILE", CELL_BC_FILE.getAbsoluteFile());
		else
			setDgeHeaderLibraryField(lib, "CELL_BC_FILE", null);
        // TODO: modify this properly to strand strategy?
        // setDgeHeaderLibraryField(lib, "USE_STRAND_INFO", USE_STRAND_INFO);
        setDgeHeaderLibraryField(lib, "RARE_UMI_FILTER_THRESHOLD", RARE_UMI_FILTER_THRESHOLD);
        setDgeHeaderLibraryField(lib, "STRAND_STRATEGY", (this.STRAND_STRATEGY == null? null: this.STRAND_STRATEGY.name()));
        setDgeHeaderLibraryField(lib, "LOCUS_FUNCTION_LIST", this.LOCUS_FUNCTION_LIST.toString());
        header.addLibrary(lib);
        header.addCommand(getCommandLine());
        final OutputStreamWriter writer = new OutputStreamWriter(out);
        new DgeHeaderCodec().encode(writer, header);
        try {
            writer.flush();
        } catch (IOException e) {
            throw new RuntimeException("Exception writing " + OUTPUT, e);
        }
    }

    private <T> void setDgeHeaderLibraryField(final DgeHeaderLibrary lib, final String key, final T value) {
        String stringValue;
        if (value != null)
			stringValue = value.toString();
		else
			stringValue = "NULL";
        lib.setTag(key, stringValue);
    }




    private void writeStats (final String gene, final Map<String, Integer> countMap, final List<String> cellBarcodes, final PrintStream out) {

        int totalCount=0;
        List<String> line = new ArrayList<>(cellBarcodes.size()+1);
        line.add(gene);
        for (String b: cellBarcodes) {
            Integer count = countMap.get(b);
            String v = "0";
            if (count!=null) {
                totalCount+=count;
                v = count.toString();
            }
            line.add(v);
        }
        if (OUTPUT_EXPRESSED_GENES_ONLY & totalCount==0) return;
        if (MIN_SUM_EXPRESSION!=null && totalCount < MIN_SUM_EXPRESSION) return;

        String h = StringUtils.join(line, "\t");
        out.println(h);
        //OutputWriterUtil.writeResult(h, out);
    }


    private void writeHeader(final PrintStream out, final List<String> cellBarcodes) {
        List<String> header = new ArrayList<>(cellBarcodes.size()+1);
        header.add("GENE");
        for (String c: cellBarcodes)
			header.add(c);
        String h = StringUtils.join(header, "\t");
        out.println(h);
    }

    public static class DESummary extends MetricBase {

        public String CELL_BARCODE;
        public int NUM_GENIC_READS;
        public int NUM_TRANSCRIPTS;
        public int NUM_GENES;

        public DESummary (final String cellBarcode) {
            this.CELL_BARCODE=cellBarcode;
            this.NUM_GENES=0;
            this.NUM_TRANSCRIPTS=0;
            this.NUM_GENIC_READS=0;
        }

    }

    static final Comparator<DESummary> TRANSCRIPT_ORDER_DESCENDING =  new Comparator<DESummary>() {
        @Override
		public int compare(final DESummary e1, final DESummary e2) {
            return (e1.NUM_TRANSCRIPTS > e2.NUM_TRANSCRIPTS ? -1 : e1.NUM_TRANSCRIPTS == e2.NUM_TRANSCRIPTS ? 0 : 1);
        }
    };

    public static Map<String, DESummary> initializeSummary(final Collection<String> cellBarcodes) {
        Map<String, DESummary> map = new HashMap<>();

        for (String s: cellBarcodes) {
            DESummary des = new DESummary(s);
            map.put(s, des);
        }
        return (map);
    }

    public static Map<String, DESummary> addToSummary(final Map<String, Integer> readCountMap, final Map<String, Integer> transcriptCountMap, final Map<String, DESummary> summaryMap) {
        for (String cellBC: transcriptCountMap.keySet()) {
            DESummary sum = summaryMap.get(cellBC);
            // for genes, it doesn't matter what the count is as long as it's > 0.  Increment by 1.
            sum.NUM_GENES++;
            // for transcripts, increment by the count.
            sum.NUM_TRANSCRIPTS+=transcriptCountMap.get(cellBC);
            sum.NUM_GENIC_READS+=readCountMap.get(cellBC);
        }
        return (summaryMap);
    }

    void writeSummary(final Collection<DESummary> summaryCollection, final File outFile) {
        MetricsFile<DESummary, Integer> out = getMetricsFile();
        List<DESummary> sc = new ArrayList<>(summaryCollection);
        Collections.sort(sc, DigitalExpression.TRANSCRIPT_ORDER_DESCENDING);
        for (DESummary z: sc)
			out.addMetric(z);
        out.write(outFile);
    }

    private SortingCollection<DGELongFormatRecord> makeSortingCollection(final Collection<String> cellBarcodeOrder) {
        final String name = getClass().getSimpleName();
        final File tmpDir = IOUtil.createTempDir(name, null);
        tmpDir.deleteOnExit();
        // 256 byte for one name, and 1/10 part of all memory for this, rough estimate
        long maxNamesInRam = Runtime.getRuntime().maxMemory() / 256 / 10;
        return SortingCollection.newInstance(
        		DGELongFormatRecord.class,
                new DGELongFormatRecordCodec(),
                new CellBarcodeOrderComparator(cellBarcodeOrder),
                (int) Math.min(maxNamesInRam, Integer.MAX_VALUE),
                tmpDir.toPath()
        );
    }


    /** Stock main method. */
    public static void main(final String[] args) {
        System.exit(new DigitalExpression().instanceMain(args));
    }

    /**
     * Collapses a bunch of strings by the edit distance.
     * If edit distance computations indicate it's greater than threshold edit distance, then the threshold is returned.
     * This is to avoid hard work on indel calculations when edit distance between two strings is high.
     * You can safely set threshold to be 3 * edit distance.
     * @param barcodes
     * @param editDistance
     * @return
     */
    /*
    public ObjectCounter <String> collapseByEditDistance (final ObjectCounter<String> barcodes, final int editDistance) {
        // map the barcode to the object so I can look up counts


        ObjectCounter <String> result = new ObjectCounter<>();
        List<String> barcodeList = barcodes.getKeysOrderedByCount(true);

        // short circuit for ED=0
        if (this.EDIT_DISTANCE==0) {
            for (String barcode: barcodeList) {
                int count=barcodes.getCountForKey(barcode);
                result.setCount(barcode, count);
            }
            return (result);
        }

        while (barcodeList.isEmpty()==false) {
            String b = barcodeList.get(0);
            barcodeList.remove(b);
            // this is still the "old" single core version.  Molecular barcode counts are small, so this may be ok.
            Set<String> closeBC = EDUtils.getInstance().getStringsWithinEditDistanceWithIndel(b,barcodeList, editDistance);
            barcodeList.removeAll(closeBC);
            // for counting.
            closeBC.add(b);
            int totalCount = 0;
            for (String bc: closeBC) {
                int count = barcodes.getCountForKey(bc);
                totalCount+=count;
            }
            result.setCount(b, totalCount);
        }
        return (result);
    }
	*/

}
