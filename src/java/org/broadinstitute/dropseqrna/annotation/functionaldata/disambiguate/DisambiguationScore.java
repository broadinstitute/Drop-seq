package org.broadinstitute.dropseqrna.annotation.functionaldata.disambiguate;

import org.broadinstitute.dropseqrna.annotation.functionaldata.FunctionalData;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import picard.annotation.LocusFunction;

import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.stream.Collectors;


public class DisambiguationScore {


    private final String cell;
    private final String molecularBarcode;

    private final ObjectCounter<String> ambiguousCount;
    private final ObjectCounter<String> antisenseCodingCount;
    private final ObjectCounter<String> senseIntronicCount;

    public DisambiguationScore(String cell, String molecularBarcode) {
        this.cell=cell;
        this.molecularBarcode=molecularBarcode;
        this.ambiguousCount=new ObjectCounter<>();
        this.antisenseCodingCount =new ObjectCounter<>();
        this.senseIntronicCount=new ObjectCounter<>();
    }

    /**
     * Track the gene(s) the read came from, and if they unambiguously assigned to an sense-intron or antisense-exon,
     * or if they assigned to both.
     * Ambiguous assignments are assigned to all gene names.
     * @param fd functional data for a single read.
     */
    public void add(List<FunctionalData> fd) {
        boolean antisenseCoding=containsAntisenseCoding(fd);
        boolean senseIntronic=containsSenseIntronic(fd);
        List<String> geneNames = getGeneNames(fd);
        // if both, then make a "metagene" of the two genes that contributed to this.
        if (antisenseCoding && senseIntronic) {
            if (geneNames.size()>2)
                throw new IllegalStateException(("More than 2 genes?"));
            for (String geneName: geneNames)
                ambiguousCount.increment(geneName);
        } else {
            if (geneNames.size()>1)
                throw new IllegalStateException(("Should only be one gene name for unambiguous reads"));
            if (antisenseCoding)
                antisenseCodingCount.increment(geneNames.getFirst());
            if (senseIntronic)
                senseIntronicCount.increment (geneNames.getFirst());
            }
    }

    /**
     * Is there evidence that at least one read for this cell/UMI is ambiguous?
     * @return True if at least one read was assigned to both sense intronic and antisense coding
     */
    public boolean hasAmbiguousReads () {
        return (this.ambiguousCount.getTotalCount() >0);
    }

    /**
     * Classify this cell/molecular barcode data as coming from one of the 4 possible categories.
     * 1. The data is unambiguous - no ambiguous reads were observed
     * 2. There was at least one ambiguous read, and no reads existed that contained clean evidence to disambiguate it
     * 3. There was at least one ambiguous read, and at least one read that could resolve the genes to either antisense
     * coding or sense intronic interpretation
     * @return
     */
    public FunctionCategory classify () {
        if (!this.hasAmbiguousReads())
            return FunctionCategory.UNAMBIGUOUS;

        // resolve by looking for gene overlaps between unambiguous and ambiguous
        if (testKeyOverlap(this.ambiguousCount, this.antisenseCodingCount))
            return FunctionCategory.ANTISENSE_CODING;

        if (testKeyOverlap(this.ambiguousCount, this.senseIntronicCount))
            return FunctionCategory.SENSE_INTRONIC;

        // if no reads resolve the data
        return FunctionCategory.AMBIGUOUS;
    }

    private boolean testKeyOverlap (ObjectCounter<String> o1, ObjectCounter<String> o2) {
        return o1.getKeys().stream().anyMatch(o2.getKeys()::contains);
    }

    private String getGeneName (List<FunctionalData> fd) {
        List<String> geneNames = getGeneNames(fd);
        String geneName=String.join(":", geneNames);
        return geneName;
    }

    private List<String> getGeneNames (List<FunctionalData> fd) {
        List<String> genes = fd.stream().map(x -> x.getGene()).collect(Collectors.toList());
        Collections.sort(genes);
        return (genes);
    }
    private boolean isAntisenseCoding (FunctionalData d) {
        return (d.getLocusFunction()== LocusFunction.CODING |d.getLocusFunction()==LocusFunction.UTR) & !d.getGeneStrand().equals(d.getReadStrand());
    }

    private boolean containsAntisenseCoding (List<FunctionalData> fd) {
        return fd.stream().anyMatch(this::isAntisenseCoding);
    }

    private boolean containsSenseIntronic (List<FunctionalData> fd) {
        return fd.stream().anyMatch(this::isSenseIntronic);
    }

    private boolean isSenseIntronic (FunctionalData d) {
        return d.getLocusFunction()==LocusFunction.INTRONIC & Objects.equals(d.getGeneStrand(), d.getReadStrand());
    }

    public String getCell() {
        return cell;
    }

    public String getMolecularBarcode() {
        return molecularBarcode;
    }

    public ObjectCounter<String> getAmbiguousCount() {
        return ambiguousCount;
    }

    public ObjectCounter<String> getAntisenseCodingCount() {
        return antisenseCodingCount;
    }

    public ObjectCounter<String> getSenseIntronicCount() {
        return senseIntronicCount;
    }

    /**
     * Get the total number of recorded reads.
     * @return
     */
    public int getTotalCount () {
        return this.ambiguousCount.getTotalCount()+this.antisenseCodingCount.getTotalCount()+this.senseIntronicCount.getTotalCount();
    }

    public String toString () {
        StringBuilder b = new StringBuilder();
        b.append("Cell [").append(this.getCell()).append("] UMI [").append(this.molecularBarcode)
                .append("] ambiguous [").append(this.ambiguousCount.toString()).append("]")
                .append("] antisense coding  [").append(this.antisenseCodingCount.toString()).append("]")
                .append("] sense intronic  [").append(this.senseIntronicCount.toString()).append("]");

        return b.toString();
    }



}
