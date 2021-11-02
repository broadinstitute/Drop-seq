package org.broadinstitute.dropseqrna.sbarro;

import htsjdk.samtools.SAMRecord;

import java.util.Map;
import java.util.Set;
import java.util.function.Predicate;

@SuppressWarnings("ClassCanBeRecord")
public class ExistingCbcUmiPredicate implements Predicate<SAMRecord> {
    private final Map<String, Set<String>> cellBarcodeUmis;

    private static final String CELL_BARCODE_TAG = "XC";
    private static final String MOLECULAR_BARCODE_TAG = "XM";

    public ExistingCbcUmiPredicate(final Map<String, Set<String>> cellBarcodeUmis) {
        this.cellBarcodeUmis = cellBarcodeUmis;
    }

    @Override
    public boolean test(SAMRecord samRecord) {
        final String cbc = samRecord.getStringAttribute(CELL_BARCODE_TAG);
        final Set<String> umis = cellBarcodeUmis.get(cbc);
        return umis != null && umis.contains(samRecord.getStringAttribute(MOLECULAR_BARCODE_TAG));
    }
}
