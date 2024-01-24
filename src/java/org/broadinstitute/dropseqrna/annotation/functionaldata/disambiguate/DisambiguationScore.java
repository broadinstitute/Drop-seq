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
    private final ObjectCounter<String> senseCodingCount;

    private boolean isComplex;

    public DisambiguationScore(String cell, String molecularBarcode) {
        this.cell=cell;
        this.molecularBarcode=molecularBarcode;
        this.ambiguousCount=new ObjectCounter<>();
        this.antisenseCodingCount =new ObjectCounter<>();
        this.senseIntronicCount=new ObjectCounter<>();
        this.senseCodingCount =new ObjectCounter<>();
        this.isComplex=false;
    }

    /**
     * Track the gene(s) the read came from, and if they unambiguously assigned to an sense-intron or antisense-exon,
     * or if they assigned to both.
     * Ambiguous assignments are assigned to all gene names.
     * @param fd functional data for a single read.
     */
    public void add(List<FunctionalData> fd) {
        List<String> geneNames = getGeneNames(fd);

        // if the read is a sense coding read that trumps all other priority, so don't evaluate
        boolean isSenseCoding = containsCoding(fd, false);
        if (isSenseCoding) {
            geneNames.forEach(senseCodingCount::increment);
            return;
        }

        boolean antisenseCoding=containsCoding(fd, true);
        boolean senseIntronic=containsSenseIntronic(fd);
        // if both, then make a "metagene" of the two genes that contributed to this.
        if (antisenseCoding && senseIntronic) {
            if (geneNames.size()>2)
                this.isComplex=true;
            geneNames.forEach(ambiguousCount::increment);
        } else {
            if (geneNames.size()!=1)
                this.isComplex=true;
            if (antisenseCoding)
                geneNames.forEach(antisenseCodingCount::increment);
            if (senseIntronic)
                geneNames.forEach(senseIntronicCount::increment);
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

        boolean overlapAntisenseCoding=testKeyOverlap(this.ambiguousCount, this.antisenseCodingCount);
        boolean overlapSenseIntronic=testKeyOverlap(this.ambiguousCount, this.senseIntronicCount);

        // if there are reads that disambiguate in both directions, it's still ambiguous!
        if (overlapAntisenseCoding && overlapSenseIntronic)
            return FunctionCategory.AMBIGUOUS;

        // resolve by looking for gene overlaps between unambiguous and ambiguous
        if (overlapAntisenseCoding)
            return FunctionCategory.ANTISENSE_CODING;

        if (overlapSenseIntronic)
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

    private boolean isCoding(FunctionalData d, boolean antisense) {
        return (d.getLocusFunction() == LocusFunction.CODING || d.getLocusFunction() == LocusFunction.UTR)
                && Objects.equals(d.getGeneStrand(), d.getReadStrand()) != antisense;
    }

    private boolean containsCoding(List<FunctionalData> fd, boolean antisense) {
        return fd.stream().anyMatch(data -> isCoding(data, antisense));
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

    public ObjectCounter<String> getSenseCodingCount() {
        return senseCodingCount;
    }

    public ObjectCounter<String> getSenseIntronicCount() {
        return senseIntronicCount;
    }

    public boolean isComplex() {
        return isComplex;
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
        b.append("Cell [").append(this.getCell()).append("] UMI [").append(this.molecularBarcode).append("] ")
                .append("Category [").append(classify().toString()).append("]")
                .append(" is complex [").append(Boolean.toString(this.isComplex)).append(" ] ")
                .append("] sense coding [").append(this.senseCodingCount.toString()).append("]")
                .append("] ambiguous [").append(this.ambiguousCount.toString()).append("]")
                .append("] antisense coding  [").append(this.antisenseCodingCount.toString()).append("]")
                .append("] sense intronic  [").append(this.senseIntronicCount.toString()).append("]");
        return b.toString();
    }



}
