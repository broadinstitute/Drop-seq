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

    // private final ObjectCounter<String> ambiguousCount;


    private final ObjectCounter<String> ambiguousAntiSenseCount;
    private final ObjectCounter<String> ambiguousSenseIntronicCount;

    private final ObjectCounter<String> antisenseCodingCount;
    private final ObjectCounter<String> senseIntronicCount;
    private final ObjectCounter<String> senseCodingCount;

    private int numReads;
    private boolean isComplex;

    public DisambiguationScore(String cell, String molecularBarcode) {
        this.cell=cell;
        this.molecularBarcode=molecularBarcode;
        // this.ambiguousCount=new ObjectCounter<>();
        this.antisenseCodingCount =new ObjectCounter<>();
        this.senseIntronicCount=new ObjectCounter<>();
        this.senseCodingCount =new ObjectCounter<>();
        this.ambiguousAntiSenseCount = new ObjectCounter<>();
        this.ambiguousSenseIntronicCount = new ObjectCounter<>();
        this.isComplex=false;
        this.numReads=0;
    }

    /**
     * Track the gene(s) the read came from, and if they unambiguously assigned to an sense-intron or antisense-exon,
     * or if they assigned to both.
     * Ambiguous assignments are assigned to all gene names.
     * @param fd functional data for a single read.
     */
    public void add(List<FunctionalData> fd) {
        this.numReads++;
        List<String> geneNames = getGeneNames(fd);

        // if the read is a sense coding read that trumps all other priority, so don't evaluate
        boolean isSenseCoding = containsCoding(fd, false);
        if (isSenseCoding) {
            geneNames.forEach(senseCodingCount::increment);
            return;
        }

        boolean antisenseCoding=containsCoding(fd, true);
        boolean senseIntronic=containsSenseIntronic(fd);
        // if both, this is the ambiguous category.
        if (antisenseCoding && senseIntronic) {
            if (geneNames.size()>2)
                this.isComplex=true;
            // more explicitly capture the ambiguous antisense coding and ambiguous sense intronic genes!
            for (FunctionalData d: fd) {
                if (isCoding(d, true))
                    this.ambiguousAntiSenseCount.increment(d.getGene());
                if (isSenseIntronic(d))
                    this.ambiguousSenseIntronicCount.increment(d.getGene());
            }
        } else {
            if (geneNames.size() != 1)
                this.isComplex = true;
            // more explicitly capture the unambiguous coding and sense intronic genes!
            for (FunctionalData d : fd) {
                if (isCoding(d, true))
                    this.antisenseCodingCount.increment(d.getGene());
                if (isSenseIntronic(d))
                    this.senseIntronicCount.increment(d.getGene());
            }
        }
    }

    /**
     * Is there evidence that at least one read for this cell/UMI is ambiguous?
     * @return True if at least one read was assigned to both sense intronic and antisense coding
     */
    public boolean hasAmbiguousReads () {
        return (this.ambiguousAntiSenseCount.getTotalCount() >0 ||  this.ambiguousSenseIntronicCount.getTotalCount() >0);
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
        if (!this.hasAmbiguousReads()) {
            if (this.senseCodingCount.getTotalCount()>0)
                return FunctionCategory.SENSE_CODING;
            if (this.antisenseCodingCount.getTotalCount()>0)
                return FunctionCategory.ANTISENSE_CODING;
            if (this.senseIntronicCount.getTotalCount()>0)
                return FunctionCategory.SENSE_INTRONIC;
        }

        boolean overlapAntisenseCoding=testKeyOverlap(this.ambiguousAntiSenseCount, this.antisenseCodingCount);
        boolean overlapSenseIntronic=testKeyOverlap(this.ambiguousSenseIntronicCount, this.senseIntronicCount);

        // if there are reads that disambiguate in both directions, it's unresolved.
        if (overlapAntisenseCoding && overlapSenseIntronic)
            return FunctionCategory.UNRESOLVED;

        // if there are no reads to disambiguate in either direction, it's unresolved
        if (this.hasAmbiguousReads() & this.senseIntronicCount.getSize()==0 & this.antisenseCodingCount.getSize()==0)
            return FunctionCategory.UNRESOLVED;

        // resolve by looking for gene overlaps between unambiguous and ambiguous
        if (overlapAntisenseCoding)
            return FunctionCategory.RESOLVED_ANTISENSE_CODING;

        if (overlapSenseIntronic)
            return FunctionCategory.RESOLVED_SENSE_INTRONIC;

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
        return this.numReads;
    }

    public ObjectCounter<String> getAmbiguousAntiSenseCount() {
        return ambiguousAntiSenseCount;
    }

    public ObjectCounter<String> getAmbiguousSenseIntronicCount() {
        return ambiguousSenseIntronicCount;
    }


    public String toString () {
        StringBuilder b = new StringBuilder();
        b.append("Cell [").append(this.getCell()).append("] UMI [").append(this.molecularBarcode).append("] ")
                .append("Category [").append(classify().toString()).append("]")
                .append(" is complex [").append(Boolean.toString(this.isComplex)).append(" ] ")
                .append("] sense coding [").append(this.senseCodingCount.toString()).append("]")
                .append("] ambiguous antisense [").append(this.ambiguousAntiSenseCount.toString()).append("]")
                .append("] ambiguous sense intronic [").append(this.ambiguousSenseIntronicCount.toString()).append("]")
                .append("] antisense coding  [").append(this.antisenseCodingCount.toString()).append("]")
                .append("] sense intronic  [").append(this.senseIntronicCount.toString()).append("]");
        return b.toString();
    }



}
