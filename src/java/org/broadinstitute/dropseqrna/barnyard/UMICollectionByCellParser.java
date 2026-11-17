package org.broadinstitute.dropseqrna.barnyard;

import htsjdk.samtools.util.*;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.utils.OrderAssertingIterator;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.File;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;

/**
 * Given a molecular barcodes file
 * (tab-separated compressed text file with header: CELL_BARCODE, GENE, MOLECULAR_BARCODE,  NUM_OBS)
 * Will yield a list of UMICollection objects (One per cell+gene pair) for each cell.
 * @author dmeyer
 */
public class UMICollectionByCellParser extends IterableOnceIterator<List<UMICollection>> {

    private final PeekableIterator<TabbedTextFileWithHeaderParser.Row> parserIter;
    private final Log log = Log.getInstance(UMICollectionByCellParser.class);
    private final ProgressLogger progress;

    private final boolean skipChimeric;

    /**
     * Constructor to parse a given molecular barcode file.
     * @param input molecular barcode file object (ends in .molBC.txt.gz)
     */
    public UMICollectionByCellParser(File input, boolean skipChimeric) {
        this(input, skipChimeric, null);
    }

    /**
     * Constructor to parse a given molecular barcode file, keeping only some of the cell barcodes.
     * @param input molecular barcode file object (ends in .molBC.txt.gz)
     * @param cellBarcodes the cell barcodes to keep; rows of any other cell barcode are skipped, and progress is logged
     *                     as a fraction of this list.  If null, all cell barcodes are kept.
     */
    public UMICollectionByCellParser(File input, boolean skipChimeric, Set<String> cellBarcodes) {
        this.skipChimeric = skipChimeric;
        this.progress = cellBarcodes == null
                ? new CellProgressLogger(log, 10000, "cell barcodes")
                : new CellProgressLogger(log, 1000, "of " + NumberFormat.getIntegerInstance().format(cellBarcodes.size()) + " cell barcodes");
        IOUtil.assertFileIsReadable(input);
        IterableOnceIterator<TabbedTextFileWithHeaderParser.Row> iterator = new OrderAssertingIterator<>(
                (new TabbedTextFileWithHeaderParser(input)).iterator(),
                Comparator.comparing((TabbedTextFileWithHeaderParser.Row o) -> o.getField(GatherMolecularBarcodeDistributionByGene.CELL_BARCODE_COLUMN))
                        .thenComparing(o -> o.getField(GatherMolecularBarcodeDistributionByGene.GENE_COLUMN))
        );
        if (cellBarcodes != null) {
            iterator = new CellBarcodeFilteringIterator(iterator, cellBarcodes);
        }
        if (skipChimeric) {
            iterator = new ChimericSkippingIterator(iterator);
        }
        parserIter = new PeekableIterator<>(iterator);
    }

    public UMICollectionByCellParser(File input) {
        this(input, false);
    }

    @Override
    public boolean hasNext() { return parserIter.hasNext(); }

    public TabbedTextFileWithHeaderParser.Row peek() { return parserIter.peek(); }

    @Override
    public List<UMICollection> next() {
        ArrayList<UMICollection> res = new ArrayList<>();

        TabbedTextFileWithHeaderParser.Row row = parserIter.next();
        UMICollection currentUMI = new UMICollection(row.getField(GatherMolecularBarcodeDistributionByGene.CELL_BARCODE_COLUMN), row.getField(GatherMolecularBarcodeDistributionByGene.GENE_COLUMN));
        currentUMI.incrementMolecularBarcodeCount(row.getField(GatherMolecularBarcodeDistributionByGene.MOLECULAR_BARCODE_COLUMN),
                Integer.parseInt(row.getField(GatherMolecularBarcodeDistributionByGene.NUM_OBS_COLUMN)));

        //when a new gene/cell is seen, make a new object and put records into that, and store the old gene/cell.
        while(parserIter.hasNext()) {
            row = parserIter.peek();
            String cell = row.getField(GatherMolecularBarcodeDistributionByGene.CELL_BARCODE_COLUMN);
            String gene = row.getField(GatherMolecularBarcodeDistributionByGene.GENE_COLUMN);
            String molBc = row.getField(GatherMolecularBarcodeDistributionByGene.MOLECULAR_BARCODE_COLUMN);

            if (!cell.equals(currentUMI.getCellBarcode())) {
                res.add(currentUMI);
                progress.record("",0);
                return res;
            }
            if (!gene.equals(currentUMI.getGeneName())) {
                res.add(currentUMI);
                currentUMI = new UMICollection(cell, gene);
            }
            currentUMI.incrementMolecularBarcodeCount(molBc, Integer.parseInt(row.getField(GatherMolecularBarcodeDistributionByGene.NUM_OBS_COLUMN)));

            parserIter.next();
        }
        res.add(currentUMI);
        if (!this.hasNext())
            progress.record("",0);
        return res;
    }

    /** Passes only the rows of the given cell barcodes. */
    private static class CellBarcodeFilteringIterator
    extends IterableOnceIterator<TabbedTextFileWithHeaderParser.Row> {
        private final Iterator<TabbedTextFileWithHeaderParser.Row> underlyingIterator;
        private final Set<String> cellBarcodes;
        private TabbedTextFileWithHeaderParser.Row nextRow = null;

        CellBarcodeFilteringIterator(final Iterator<TabbedTextFileWithHeaderParser.Row> iterator, final Set<String> cellBarcodes) {
            this.underlyingIterator = iterator;
            this.cellBarcodes = cellBarcodes;
        }

        @Override
        public boolean hasNext() {
            while (nextRow == null && underlyingIterator.hasNext()) {
                final TabbedTextFileWithHeaderParser.Row row = underlyingIterator.next();
                if (cellBarcodes.contains(row.getField(GatherMolecularBarcodeDistributionByGene.CELL_BARCODE_COLUMN))) {
                    nextRow = row;
                }
            }
            return nextRow != null;
        }

        @Override
        public TabbedTextFileWithHeaderParser.Row next() {
            if (!hasNext()) {
                throw new NoSuchElementException();
            }
            final TabbedTextFileWithHeaderParser.Row row = nextRow;
            nextRow = null;
            return row;
        }
    }

    private static class ChimericSkippingIterator
    extends IterableOnceIterator<TabbedTextFileWithHeaderParser.Row> {
        private final PeekableIterator<TabbedTextFileWithHeaderParser.Row> underlyingIterator;

        public ChimericSkippingIterator(Iterator<TabbedTextFileWithHeaderParser.Row> iterator) {
            underlyingIterator = new PeekableIterator<>(iterator);
        }

        private void skipChimeric() {
            while (underlyingIterator.hasNext() && underlyingIterator.peek().getField("CHIMERIC").equals("true")) {
                underlyingIterator.next();
            }
        }

        @Override
        public boolean hasNext() {
            skipChimeric();
            return underlyingIterator.hasNext();
        }

        @Override
        public TabbedTextFileWithHeaderParser.Row next() {
            skipChimeric();
            return underlyingIterator.next();
        }

    }

    /**
     * Logs progress every 10,000 cells.  ProgressLogger was written for records with a genomic position and ends its message
     * with "Last read position: chrom:pos", which is meaningless here (there is no position in a text report), so that is left off.
     */
    private static class CellProgressLogger extends ProgressLogger {
        private final Log log;

        CellProgressLogger(final Log log, final int interval, final String noun) {
            super(log, interval, "Processed", noun);
            this.log = log;
        }

        @Override
        protected void log(final String... message) {
            log.info(String.join("", message).replaceFirst("\\s*Last read position:.*$", ""));
        }
    }
}

