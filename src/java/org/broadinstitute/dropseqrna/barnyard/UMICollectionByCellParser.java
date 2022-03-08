package org.broadinstitute.dropseqrna.barnyard;

import htsjdk.samtools.util.*;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.utils.OrderAssertingIterator;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

/**
 * Given a molecular barcodes file
 * (tab-separated compressed text file with header: Cell Barcode, Gene, Molecular_Barcode,  Num_Obs)
 * Will yield a list of UMICollection objects (One per cell+gene pair) for each cell.
 * @author dmeyer
 */
public class UMICollectionByCellParser extends IterableOnceIterator<List<UMICollection>> {

    private PeekableIterator<TabbedTextFileWithHeaderParser.Row> parserIter;
    private Log log = Log.getInstance(UMICollectionByCellParser.class);
    private ProgressLogger progress = new ProgressLogger(log, 100, "Processed", "cells");

    /**
     * Constructor to parse a given molecular barcode file.
     * @param input molecular barcode file object (ends in .molBC.txt.gz)
     */
    UMICollectionByCellParser(File input) {
        IOUtil.assertFileIsReadable(input);
        parserIter = new PeekableIterator<>(
                new OrderAssertingIterator<>(
                        (new TabbedTextFileWithHeaderParser(input)).iterator(),
                        Comparator.comparing((TabbedTextFileWithHeaderParser.Row o) -> o.getField("Cell Barcode"))
                                .thenComparing(o -> o.getField("Gene"))
                )
        );
    }

    @Override
    public boolean hasNext() { return parserIter.hasNext(); }

    public TabbedTextFileWithHeaderParser.Row peek() { return parserIter.peek(); }

    @Override
    public List<UMICollection> next() {
        ArrayList<UMICollection> res = new ArrayList<>();

        TabbedTextFileWithHeaderParser.Row row = parserIter.next();
        UMICollection currentUMI = new UMICollection(row.getField("Cell Barcode"), row.getField("Gene"));
        currentUMI.incrementMolecularBarcodeCount(row.getField("Molecular_Barcode"),
                Integer.parseInt(row.getField("Num_Obs")));

        //when a new gene/cell is seen, make a new object and put records into that, and store the old gene/cell.
        while(parserIter.hasNext()) {
            row = parserIter.peek();
            String cell = row.getField("Cell Barcode");
            String gene = row.getField("Gene");
            String molBc = row.getField("Molecular_Barcode");

            if (!cell.equals(currentUMI.getCellBarcode())) {
                res.add(currentUMI);
                progress.record("",0);
                return res;
            }
            if (!gene.equals(currentUMI.getGeneName())) {
                res.add(currentUMI);
                currentUMI = new UMICollection(cell, molBc);
            }
            currentUMI.incrementMolecularBarcodeCount(molBc, Integer.parseInt(row.getField("Num_Obs")));

            parserIter.next();
        }
        res.add(currentUMI);
        if (!this.hasNext())
            progress.record("",0);
        return res;
    }
}

