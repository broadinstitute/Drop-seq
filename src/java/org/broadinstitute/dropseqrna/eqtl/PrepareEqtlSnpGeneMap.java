/*
 * MIT License
 *
 * Copyright 2020 Broad Institute
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
package org.broadinstitute.dropseqrna.eqtl;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

@CommandLineProgramProperties(
        summary = "Generate list of {SNP, gene} pairs to be considered in eQTL discovery.  " +
                "There are 3 methods for identifying pairs: 1) SNP is within CIS_DIST of gene; " +
                "2) SNP is in an interval in NS_INTERVAL_LIST, and is within NS_CIS_DIST of gene;" +
                "3) SNP is in a gene-specific interval in GS_INTERVAL_LIST, and is within GS_CIST_DIST of that gene.  " +
                "At least one method must be specified.  If a {SNP, gene} pair satisfied any method" +
                "specified, it will be written to the output.",
        oneLineSummary = "Generate list of {SNP, gene} pairs to be considered in eQTL discovery.  ",
        programGroup = DropSeq.class
)
public class PrepareEqtlSnpGeneMap  extends CommandLineProgram {
    private final String FIELD_SEPARATOR = "\t";
    private final String SNP_HEADER = "snp";
    private final String GENE_HEADER = "gene";
    private final String SOURCE_HEADER = "source";
    private final Log LOG = Log.getInstance(PrepareEqtlSnpGeneMap.class);

    @Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
    doc="Map of SNPs and genes to be considered together in eQTL discovery")
    public File OUTPUT;

    @Argument(doc="The SNPs to be considered.")
    public File SNP_LOCATIONS;

    @Argument (doc="The genes to be considered.")
    public File GENE_LOCATION_FILE;

    @Argument(doc="Consider a SNP this close to a gene.  If not set, do not consider SNPs by distance only.",
            optional = true)
    public Integer CIS_DIST;

    @Argument(doc="List of non-gene-specific intervals to find SNPs to be tested with all genes within NS_CIS_DIST of the SNP.",
            optional = true)
    public File NS_INTERVAL_LIST;

    @Argument(doc="Test genes this close to SNPs in intervals in NS_INTERVAL_LIST.")
    public int NS_CIS_DIST = Integer.MAX_VALUE;

    @Argument(doc="List of gene-specific intervals to find SNPs to be tested with the gene within GS_CIS_DIST of the SNP. " +
            "The name of the interval is the gene for which SNPs in the interval are to be tested.",
            optional = true)
    public File GS_INTERVAL_LIST;

    @Argument(doc="Test genes this close to SNPs in intervals in GS_INTERVAL_LIST.")
    public int GS_CIS_DIST = Integer.MAX_VALUE;

    // When invoked from PrepareEqtlData, getCommandLine() returns null, so this holds the
    // PrepareEqtlData command line
    private String commandLine;

    public static void main(final String[] args) {
        new PrepareEqtlSnpGeneMap().instanceMainWithExit(args);
    }

    public PrepareEqtlSnpGeneMap() {
    }

    public PrepareEqtlSnpGeneMap(File OUTPUT, File SNP_LOCATIONS, File GENE_LOCATION_FILE, Integer CIS_DIST,
                                 File NS_INTERVAL_LIST, int NS_CIS_DIST, File GS_INTERVAL_LIST, int GS_CIS_DIST,
                                 String commandLine) {
        this.OUTPUT = OUTPUT;
        this.SNP_LOCATIONS = SNP_LOCATIONS;
        this.GENE_LOCATION_FILE = GENE_LOCATION_FILE;
        this.CIS_DIST = CIS_DIST;
        this.NS_INTERVAL_LIST = NS_INTERVAL_LIST;
        this.NS_CIS_DIST = NS_CIS_DIST;
        this.GS_INTERVAL_LIST = GS_INTERVAL_LIST;
        this.GS_CIS_DIST = GS_CIS_DIST;
        this.commandLine = commandLine;
    }

    void processData() {
        final String[] validationErrors = customCommandLineValidation();
        if (validationErrors != null) {
            throw new RuntimeException("Command line validation errors: " +
                    StringUtil.join("\n", validationErrors));
        }
        if (doWork() != 0) {
            throw new RuntimeException("doWork returned non-zero");
        }
    }

    @Override
    protected String[] customCommandLineValidation() {
        final ArrayList<String> list = new ArrayList<>();
        if (this.CIS_DIST != null && this.CIS_DIST < 0) {
            list.add("CIS_DIST cannot be negative.");
        }
        if (this.NS_CIS_DIST < 0) {
            list.add("NS_CIS_DIST cannot be negative.");
        }
        if (this.GS_CIS_DIST < 0) {
            list.add("GS_CIS_DIST cannot be negative.");
        }
        if (this.CIS_DIST == null && this.NS_INTERVAL_LIST == null && this.GS_INTERVAL_LIST == null) {
            list.add("At least one of CIS_DIST, NS_INTERVAL_LIST, GS_INTERVAL_LIST must be specified");
        }
        return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);
    }

    @Override
    protected int doWork() {
        if (commandLine == null) {
            commandLine = getCommandLine();
        }
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsReadable(SNP_LOCATIONS);
        IOUtil.assertFileIsReadable(GENE_LOCATION_FILE);
        if (NS_INTERVAL_LIST!=null) IOUtil.assertFileIsReadable(NS_INTERVAL_LIST);
        if (GS_INTERVAL_LIST!=null) IOUtil.assertFileIsReadable(GS_INTERVAL_LIST);

        final List<Interval> snps = loadSNPs();
        final List<Interval> genes = loadGenes();

        // Cap the various cis dists at a conservative maximum, to avoid integer overflow.
        // Could do this by chromosome for more accuracy, but this should be fine.
        final int smallestSnpPos = snps.stream().mapToInt(Interval::getStart).min().orElse(0);
        final int largestSnpPos = snps.stream().mapToInt(Interval::getEnd).max().orElse(0);
        final int smallestGenePos = genes.stream().mapToInt(Interval::getStart).min().orElse(0);
        final int largestGenePos = genes.stream().mapToInt(Interval::getEnd).max().orElse(0);
        final int maxCisDist = Math.max(largestGenePos - smallestSnpPos, largestSnpPos - smallestGenePos);
        final int cisDist = (CIS_DIST != null? Math.min(maxCisDist, CIS_DIST): -1);
        final int nsCisDist = Math.min(maxCisDist, NS_CIS_DIST);
        final int gsCisDist = Math.min(maxCisDist, GS_CIS_DIST);

        final List<CisDistGeneFinder> geneFinders = new ArrayList<>(3);
        if (CIS_DIST != null) {
            geneFinders.add(new CisDistGeneFinder(genes, cisDist));
        }
        if (NS_INTERVAL_LIST != null) {
            geneFinders.add(new NonSpecificIntervalGeneFinder(genes, nsCisDist, NS_INTERVAL_LIST));
        }
        if (GS_INTERVAL_LIST != null) {
            geneFinders.add(new GeneSpecificIntervalGeneFinder(genes, gsCisDist, GS_INTERVAL_LIST));
        }

        // For each SNP in SNP_LOCATIONS, get the set of genes to map to using all enabled methods.
        PrintStream snpGeneMapWriter = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(OUTPUT));
        snpGeneMapWriter.println("# " + commandLine);
        snpGeneMapWriter.println(StringUtil.join(FIELD_SEPARATOR, SNP_HEADER, GENE_HEADER, SOURCE_HEADER));
        final ProgressLogger progressLogger = new ProgressLogger(LOG, 1000000, "processed", "SNPs");
        for (final Interval snp: snps) {
            progressLogger.record(snp.getContig(), snp.getStart());
            final Multimap<String, String> genesForSnp = HashMultimap.create();
            for (final CisDistGeneFinder geneFinder: geneFinders) {
                Set<String> theseGenes = geneFinder.getCandidateGenesForSnp(snp);
                if (!theseGenes.isEmpty()) {
                    LOG.debug(String.format("%s: %s %s",
                            geneFinder.getClass().getSimpleName(), snp.getName(), StringUtil.join(",", theseGenes)));
                }
                for (final String gene: theseGenes) {
                    genesForSnp.put(gene, geneFinder.sourceLabel());
                }
            }
            for (final String gene: genesForSnp.keySet()) {
                final String sources = StringUtil.join(",", genesForSnp.get(gene));
                snpGeneMapWriter.println(StringUtil.join(FIELD_SEPARATOR, snp.getName(), gene, sources));
            }
        }
        snpGeneMapWriter.close();
        LOG.info("Done writing SNP-gene map");
        return 0;
    }

    private List<Interval> loadSNPs() {
        final TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(SNP_LOCATIONS);
        final ArrayList<Interval> ret = new ArrayList<>();
        for (final TabbedTextFileWithHeaderParser.Row row: parser) {
            final String snp = row.getField("snp");
            final String chr = row.getField("chr");
            final int pos = row.getIntegerField("pos");
            final int end = row.getIntegerField("end");
            ret.add(new Interval(chr, pos, end, false, snp));
        }
        CloserUtil.close(parser);
        if (ret.isEmpty()) {
            throw new RuntimeException(SNP_LOCATIONS.getAbsolutePath() + " is empty");
        }
        return ret;
    }

    private List<Interval> loadGenes() {
        final TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(GENE_LOCATION_FILE);
        final ArrayList<Interval> ret = new ArrayList<>();
        for (final TabbedTextFileWithHeaderParser.Row row: parser) {
            final String gene = row.getField("geneid");
            final String chr = row.getField("chr");
            final int start = row.getIntegerField("s1");
            final int end = row.getIntegerField("s2");
            ret.add(new Interval(chr, start, end, false, gene));
        }
        CloserUtil.close(parser);
        if (ret.isEmpty()) {
            throw new RuntimeException(GENE_LOCATION_FILE.getAbsolutePath() + " is empty");
        }
        return ret;
    }

    /**
     * Base class that finds genes that are within cisDist of given SNP
     */
    private class CisDistGeneFinder {
        private final OverlapDetector<String> geneOverlapDetector;

        public CisDistGeneFinder(final Iterable<Interval> genes, final int cisDist) {
            geneOverlapDetector = new OverlapDetector<>(-cisDist, 0);
            for (final Interval gene: genes) {
                geneOverlapDetector.addLhs(gene.getName(), gene);
            }
        }

        public Set<String> getCandidateGenesForSnp(Interval snp) {
            return geneOverlapDetector.getOverlaps(snp);
        }

        public String sourceLabel() {
            return "CD";
        }
    }

    /**
     * Tests if SNP is in one of the intervals given, and if so returns genes with cisDist of SNP.
     */
    private class NonSpecificIntervalGeneFinder
            extends CisDistGeneFinder {
        final OverlapDetector<Interval> overlapDetector;

        public NonSpecificIntervalGeneFinder(final Iterable<Interval> genes, final int cisDist,
                                             final File nonSpecificIntervalList) {
            super(genes, cisDist);
            overlapDetector = OverlapDetector.create(IntervalList.fromFile(nonSpecificIntervalList).getIntervals());
        }

        @Override
        public String sourceLabel() {
            return "NS";
        }

        @Override
        public Set<String> getCandidateGenesForSnp(Interval snp) {
            if (overlapDetector.overlapsAny(snp)) {
                return super.getCandidateGenesForSnp(snp);
            } else {
                return Collections.emptySet();
            }
        }
    }

    /**
     * Searches gene-specific interval list for candidate genes, then checks that SNP is within
     * cisDist of gene.
     */
    private class GeneSpecificIntervalGeneFinder
            extends CisDistGeneFinder {
        final OverlapDetector<String> overlapDetector;

        public GeneSpecificIntervalGeneFinder(final Iterable<Interval> genes, final int cisDist,
                                             final File nonSpecificIntervalList) {
            super(genes, cisDist);
            overlapDetector = new OverlapDetector<>(0, 0);
            for (final Interval interval : IntervalList.fromFile(nonSpecificIntervalList).getIntervals()) {
                overlapDetector.addLhs(interval.getName(), interval);
            }
        }

        @Override
        public String sourceLabel() {
            return "GS";
        }

        @Override
        public Set<String> getCandidateGenesForSnp(Interval snp) {
            final Set<String> ret = new HashSet<>(super.getCandidateGenesForSnp(snp));
            ret.retainAll(overlapDetector.getOverlaps(snp));
            return ret;
        }

    }
}
