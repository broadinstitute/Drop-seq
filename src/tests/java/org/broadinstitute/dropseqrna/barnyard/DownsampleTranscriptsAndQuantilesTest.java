package org.broadinstitute.dropseqrna.barnyard;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.BarcodeSimulator;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.utils.TestUtils;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.testng.Assert;
import org.testng.annotations.Test;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import static com.google.common.collect.Streams.forEachPair;

/**
 * DownsampleTranscriptsAndQuantilesTest
 * Testing suite for DownsampleTranscriptsAndQuantiles
 */
public class DownsampleTranscriptsAndQuantilesTest {

    private static final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/transcriptome/barnyard");
    private final File umiCollectionFile = new File(TEST_DATA_DIR,"UMICollectionFile.txt.gz");

    private UMICollection makeUMICollection(String cell, String gene, String[] mbcounts) {
        // Shorthand for making UMICollection objects
        UMICollection res = new UMICollection(cell, gene);
        for (String s : mbcounts) {
            String[] xx = s.split("=");
            res.incrementMolecularBarcodeCount(xx[0], Integer.parseInt(xx[1]));
        }
        return res;
    }

    @Test
    public void testUMICollectionParser() {
        // Read in test file and make sure that we can parse things alright. Manually validate line by line.
        UMICollectionByCellParser iter = new UMICollectionByCellParser(umiCollectionFile);

        // Manually created UMICollection objects to validate umiCollectionFile
        String[] cellBcs= {"AACCTATGGGCC","AACCTCGACGCC","ACGGGCCGAACG","ACTGCCTTGGAA","AGACTTCCCTCG","AGGGAAAATTGA",
                "ATCAGGGACAGA","GCGCAGAGATAC","TCGACCTGTACT","TGGCGAAGAGAT","TTATTATGGCCT", "TTGCCTTACGCG"};
        String[] genes = {"HUMAN_10:100175955-100206684:HPS1", "HUMAN_10:100175955-100206684:HPS1",
                "HUMAN_10:100175955-100206684:HPS1","HUMAN_10:100175955-100206684:HPS1",
                "HUMAN_10:100175955-100206684:HPS1","HUMAN_10:100175955-100206684:HPS1",
                "HUMAN_10:100143322-100174941:PYROXD2","HUMAN_10:100143322-100174941:PYROXD2",
                "HUMAN_10:100143322-100174941:PYROXD2","HUMAN_10:100143322-100174941:PYROXD2",
                "HUMAN_10:100007447-100028007:LOXL4","HUMAN_10:100143322-100174941:PYROXD2"};
        String[] mbcounts = {"TTGTTGTG=22,TGGGGCAG=39,GTGTGGAT=1", "ACTGCCTC=12,GTGAGGCC=22",
                "ACCTGGGC=15,TGGGGAGG=17,GGGAAGTT=46,GTGGTGCA=1,GGCAGGTG=7,CAGTAGGC=33","GACTGGGA=17","GCATCACG=2",
                "TGGCTTAT=26","GGAGGTTT=24","TATGGTCC=18","GATGGGGG=54","GTTTTGGG=18","AACGGACG=6","TATGGCAT=1"};

        List<UMICollection> manualValidationLines = IntStream.range(0,cellBcs.length)
            .mapToObj(i -> makeUMICollection(cellBcs[i], genes[i], mbcounts[i].split(",")))
            .collect(Collectors.toList());

        int i = 0;
        while (iter.hasNext()) {
            List<UMICollection> a = iter.next();
            UMICollection c1 = a.get(0);
            UMICollection c2 = manualValidationLines.get(i);
            Assert.assertEquals(c1.getCellBarcode(), c2.getCellBarcode());
            Assert.assertEquals(c1.getMolecularBarcodeCounts(), c2.getMolecularBarcodeCounts());
            i++;
        }
        CloserUtil.close(iter);
    }

   private void writeSyntheticData(File file) {
        PrintStream writer = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(file));
        writer.println(String.join("\t", GatherMolecularBarcodeDistributionByGene.CELL_BARCODE_COLUMN, GatherMolecularBarcodeDistributionByGene.GENE_COLUMN, 
        		GatherMolecularBarcodeDistributionByGene.MOLECULAR_BARCODE_COLUMN, GatherMolecularBarcodeDistributionByGene.NUM_OBS_COLUMN));
        int n = 1000;
        long seed = 140L;
        Random rando = new Random(seed);
        int[] numObs = rando.ints(n, 1, 10).toArray();
        BarcodeSimulator bs = new BarcodeSimulator(14, seed);
        List<String> cellbcs = bs.getRandomBarcodes(n).stream().sorted().collect(Collectors.toList());
        List<String> molbcs =  bs.getRandomBarcodes(n);
        for (int i = 0; i < n; i++) {
            int w = rando.ints(1,1,100).toArray()[0];
            for (int j=i; j < i+w; j++) {
                if (j >= n) break;
                writer.println(String.join("\t", cellbcs.get(i), "genie", molbcs.get(j), String.valueOf(numObs[j])));
            }
            i+=w-1;
        }
        writer.flush();
        writer.close();
    }

    private List<String> getCol(File file, String columnName) {
        TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(file);
        List<String> res = new ArrayList<>();
        for (TabbedTextFileWithHeaderParser.Row row : parser) {
            res.add(row.getField(columnName));
        }
        CloserUtil.close(parser);
        return res;
    }

    private List<Integer> getIdentityCol(File file) {
        TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(file);
        List<Integer> res = new ArrayList<>();
        for (TabbedTextFileWithHeaderParser.Row row : parser) {
            res.add(row.getIntegerField("1"));
        }
        CloserUtil.close(parser);
        return res;
    }

    @Test
    public void testDoWorkSynthetic() {
        // Testing the pipes
        DownsampleTranscriptsAndQuantiles c = new DownsampleTranscriptsAndQuantiles();
        String prefix = "DownsampleTranscriptsAndQuantilesTest";
        c.RANDOM_SEED = 140;
        c.INPUT = TestUtils.getTempReportFile(prefix, ".input.txt.gz");
        writeSyntheticData(c.INPUT);
        try {
            c.OUTPUT_DOWNSAMPLING_FILE = File.createTempFile(prefix, "JavaVersion.downsampling.txt");
            c.OUTPUT_DOWNSAMPLING_FILE.deleteOnExit();
            c.OUTPUT_QUANTILE_FILE = File.createTempFile(prefix, "JavaVersion.quantiles.txt");
            c.OUTPUT_QUANTILE_FILE.deleteOnExit();
        } catch (Exception e) {
            Assert.fail();
        }

        int ret = c.doWork();
        Assert.assertEquals(ret, 0);
        Assert.assertEquals(getIdentityCol(c.OUTPUT_DOWNSAMPLING_FILE).stream().sorted().collect(Collectors.toList()),
                c.getTranscriptsPerCell().getCounts().stream().sorted().collect(Collectors.toList()));
    }

    @Test
    public void testDoWorkSyntheticMultiThreaded() {
        // Testing the pipes
        DownsampleTranscriptsAndQuantiles c = new DownsampleTranscriptsAndQuantiles();
        String prefix = "DownsampleTranscriptsAndQuantilesTest";
        c.RANDOM_SEED = 140;
        c.NUM_THREADS = 2;
        c.INPUT = TestUtils.getTempReportFile(prefix, ".input.txt.gz");
        writeSyntheticData(c.INPUT);
        try {
            c.OUTPUT_DOWNSAMPLING_FILE = File.createTempFile(prefix, "JavaVersion.downsampling.txt");
            c.OUTPUT_DOWNSAMPLING_FILE.deleteOnExit();
            c.OUTPUT_QUANTILE_FILE = File.createTempFile(prefix, "JavaVersion.quantiles.txt");
            c.OUTPUT_QUANTILE_FILE.deleteOnExit();
        } catch (Exception e) {
            Assert.fail();
        }

        int ret = c.doWork();
        Assert.assertEquals(ret, 0);
        Assert.assertEquals(getIdentityCol(c.OUTPUT_DOWNSAMPLING_FILE).stream().sorted().collect(Collectors.toList()),
                c.getTranscriptsPerCell().getCounts().stream().sorted().collect(Collectors.toList()));
    }

    private void downsampleInputFile(File infile, File outfile, double downsampleRate) {
        TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(infile);
        PrintStream writer = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(outfile));
        writer.println(String.join("\t", GatherMolecularBarcodeDistributionByGene.CELL_BARCODE_COLUMN, GatherMolecularBarcodeDistributionByGene.GENE_COLUMN, 
        		GatherMolecularBarcodeDistributionByGene.MOLECULAR_BARCODE_COLUMN, GatherMolecularBarcodeDistributionByGene.NUM_OBS_COLUMN));
        Random random = new Random(4L);
        for (TabbedTextFileWithHeaderParser.Row row : parser) {
            if (random.nextDouble() < downsampleRate) {
                writer.println(String.join("\t", 
                		row.getField(GatherMolecularBarcodeDistributionByGene.CELL_BARCODE_COLUMN), 
                		row.getField(GatherMolecularBarcodeDistributionByGene.GENE_COLUMN),
                        row.getField(GatherMolecularBarcodeDistributionByGene.MOLECULAR_BARCODE_COLUMN), 
                        row.getField(GatherMolecularBarcodeDistributionByGene.NUM_OBS_COLUMN)));
            }
        }
        writer.flush();
        writer.close();
    }

    private void runDownsampleTranscriptsAndQuantiles(File infile, File cellBcsFile,
                                                      File outfileDownsampling, File outfileQuantiles) {
        DownsampleTranscriptsAndQuantiles c = new DownsampleTranscriptsAndQuantiles();
        c.INPUT = infile;
        c.CELL_BC_FILE = cellBcsFile;
        c.RANDOM_SEED = 140;
        c.OUTPUT_DOWNSAMPLING_FILE = outfileDownsampling;
        c.OUTPUT_QUANTILE_FILE = outfileQuantiles;
        int ret = c.doWork();
        Assert.assertEquals(ret, 0);
    }

    private void integrationTestDownsampling(File javaDownsamplingFile, File rDownsamplingFile) {
        TabbedTextFileWithHeaderParser jp = new TabbedTextFileWithHeaderParser(javaDownsamplingFile);
        TabbedTextFileWithHeaderParser rp = new TabbedTextFileWithHeaderParser(rDownsamplingFile);
        List<String> columnNames = new ArrayList<>(jp.columnLabelsList());
        columnNames.remove("CELL_BARCODE");
        Assert.assertEquals(columnNames, rp.columnLabelsList());


        HashMap<String,SimpleRegression> regressions = new HashMap<>();
        columnNames.forEach(x->regressions.put(x, new SimpleRegression()));

        Iterator<TabbedTextFileWithHeaderParser.Row> rpi = rp.iterator();
        for (TabbedTextFileWithHeaderParser.Row jr : jp) {
            final TabbedTextFileWithHeaderParser.Row rr = rpi.next();
            // For each column name, add value to appropriate regression
            columnNames.forEach(d -> regressions.get(d).addData(jr.getIntegerField(d), rr.getIntegerField(d)));
        }
        for (SimpleRegression reg : regressions.values()) {
            Assert.assertTrue(reg.getRSquare() > 0.9);
        }
    }

    private void integrationTestQuantiles(File javaQuantileFile, File rQuantileFile) {
        // Make sure similar quantiles files get produced
        List<String> javaQuantileCol = getCol(javaQuantileFile, "quantile");
        List<String> rQuantileCol = getCol(rQuantileFile, "quantile");
        Assert.assertEquals(javaQuantileCol, rQuantileCol);

        Stream<Double> s1 = getCol(javaQuantileFile, "cumulative_num_cells").stream().map(Double::parseDouble);
        Stream<Double> s2 = getCol(rQuantileFile, "cumulative_num_cells").stream().map(Double::parseDouble);
        SimpleRegression regression = new SimpleRegression();
        forEachPair(s1, s2, regression::addData);
        Assert.assertTrue(regression.getRSquare() > 0.9);

        s1 = getCol(javaQuantileFile, "median_transcripts").stream().map(Double::parseDouble);
        s2 = getCol(rQuantileFile, "median_transcripts").stream().map(Double::parseDouble);
        forEachPair(s1, s2, regression::addData);
        Assert.assertTrue(regression.getRSquare() > 0.9);
    }

    @Test
    public void integrationTest() throws IOException {
        /* Integration test against output from R code to
         * To recreate R outputs to test against, run the following in R from the dropseqrna github directory:
         *
         * library(DropSeq.barnyard)
         * downsampleTranscriptsAndQuantiles(
         *     outDownsamplingFile="transcriptome_java/testdata/org/broadinstitute/transcriptome/barnyard/downsampled.d35Ngn2plusGlia_E7.Rversion.downsampling.txt",
         *     outQuantileFile="transcriptome_java/testdata/org/broadinstitute/transcriptome/barnyard/downsampled.d35Ngn2plusGlia_E7.Rversion.quantiles.txt",
         *     molecularBarcodeDistributionByGeneFile="transcriptome_java/testdata/org/broadinstitute/transcriptome/barnyard/downsampled.d35Ngn2plusGlia_E7.molBC.txt.gz",
         *     selectedCellsFile="transcriptome_java/testdata/org/broadinstitute/transcriptome/barnyard/downsampled.d35Ngn2plusGlia_E7.selectedCellBarcodes.txt",
         *     random.seed=140)
         */
        final String prefix = "downsampled.d35Ngn2plusGlia_E7";
        final File inputFile = new File(TEST_DATA_DIR,prefix+".molBC.txt.gz");
        final File cellBcsFile = new File(TEST_DATA_DIR, prefix+ ".selectedCellBarcodes.txt");
        final File rDownsamplingFile = new File(TEST_DATA_DIR, prefix+".Rversion.downsampling.txt");
        final File javaDownsamplingFile = File.createTempFile(prefix, ".JavaVersion.downsampling.txt");
        javaDownsamplingFile.deleteOnExit();
        final File javaQuantileFile = File.createTempFile(prefix, ".JavaVersion.quantiles.txt");
        javaQuantileFile.deleteOnExit();
        final File rQuantileFile = new File(TEST_DATA_DIR, prefix+".Rversion.quantiles.txt");
        runDownsampleTranscriptsAndQuantiles(inputFile, cellBcsFile, javaDownsamplingFile, javaQuantileFile);
        integrationTestDownsampling(javaDownsamplingFile, rDownsamplingFile);
        integrationTestQuantiles(javaQuantileFile, rQuantileFile);
    }
}
