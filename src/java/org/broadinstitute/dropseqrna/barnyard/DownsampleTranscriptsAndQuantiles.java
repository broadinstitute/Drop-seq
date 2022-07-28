/*
 * MIT License
 *
 * Copyright 2019 Broad Institute
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

import com.google.common.collect.Lists;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.UMICollection;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import picard.cmdline.CommandLineProgram;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.function.Function;
import java.util.function.ToIntFunction;
import java.util.stream.Collectors;

import static com.google.common.math.Quantiles.percentiles;

/** DownsampleTranscriptsAndQuantiles
 * @author dmeyer
 */
@CommandLineProgramProperties(
        summary = "For each cell, estimate the number of UMIs that we would observe if we were to sequence at "+
                "different fractions of read depth, or downsampling rates (Defaults: {.1,.2,...,.9,1}).",
        oneLineSummary = "Estimates UMIs per cell barcode if we were to have fewer reads.",
        programGroup = DropSeq.class)
public class DownsampleTranscriptsAndQuantiles extends CommandLineProgram {
    /* Inputs */
    @Argument(doc="Compressed tab-separated text file with a header with columns 'CELL_BARCODE', "+
            "'GENE', 'MOLECULAR_BARCODE', and 'NUM_OBS'. Each row represents a unique molecular barcode (UMI)," +
            "and, respectively, the columns represent the oligonucleotide cell barcode sequence associated with the " +
            "transcript, the gene to which the transcript best aligned, the oligonucleotide molecular barcode " +
            "sequence, and the number of observations of that UMI. File MUST be in ascending alphabetic order " +
            "by cell barcode then by gene.")
    public File INPUT;

    @Argument(doc="Unordered headerless text file containing newline-separated cell barcode sequences to include in " +
            "downsampling analysis. Other cell barcodes in input file will be ignored.")
    public File CELL_BC_FILE=null;

    @Argument(doc = "Random seed to use if deterministic behavior is desired. " +
            "Setting to null will cause multiple invocations to produce different results.")
    public Integer RANDOM_SEED = 1;

    @Argument(doc="List of rates (between 0.0 and 1.0) to test. For each cell, all of these rates will be " +
            "used to test how many transcripts the cell would likely retain if given these different rates " +
            "of read depth.")
    public List<Double> DOWNSAMPLING_RATES = Arrays.asList(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0);

    @Argument(doc="List of quantiles (as percentage) at which to measure median number of transcripts per cell and "+
            "cumulative number of cells that fall within or below that quantile." )
    public List<Integer> QUANTILES = Arrays.asList(0,1,10,20,30,40,50,60,70,80,90,99,100);

    /* Outputs */
    @Argument(doc="Main output of program. Tab-separated text file with a header containing downsampled transcript "+
            "counts per cell at different downsampling rates. " +
            "Rows represent integer counts of the number of transcripts assigned to the cell at different "+
            "downsampling rates (in the order of the provided DOWNSAMPLING RATES).")
    public File OUTPUT_DOWNSAMPLING_FILE;

    @Argument(doc="Tab separated file with header containing quantiles (as specified in QUANTILES) "+
            "of downsampled data. Summary of overall distribution of transcripts across cells. " +
            "Can be derived using the column of OUTPUT_DOWNSAMPLING_FILE representing a downsampling rate of 1.")
    public File OUTPUT_QUANTILE_FILE;

    @Argument(doc="Number of threads to use.  If positive, that number of threads are used.  If non-positive, " +
            "the number of processors + NUM_THREADS is used.", optional = true)
    public int NUM_THREADS=1;

    private Random random;
    private final ObjectCounter<String> transcriptsPerCell = new ObjectCounter<>();
    private Set<String> cellBarcodes;
    private ForkJoinPool pool;

    /**
     * Consumes a list of UMICollection objects, all with the same cell barcode, and returns a list of integers that
     * represents (in order) the number of transcripts (unique UMIs per cell) at each downsampling rate specified by
     * command line argument DOWNSAMPLING_RATES.
     *
     * @param umiCollections list of UMICollection objects representing transcripts of a single cell
     * @return the number of distinct molecular barcodes for a particular cell retained at different downsampling rates
     */
    private List<Integer> getDownsampledCounts(List<UMICollection> umiCollections) {
        Function<Double, ToIntFunction<UMICollection>> downsampledCountFunGen = rate -> umi -> {
            ObjectCounter<String> downsampledCounts = umi.getDownsampledMolecularBarcodeCounts(rate, random);
            downsampledCounts.filterByMinCount(1);
            return downsampledCounts.getSize();
        };
        return DOWNSAMPLING_RATES.stream()
                .map(rate -> {
                    try {
                        if (NUM_THREADS == 1) { // For reproducible behavior
                            return umiCollections.stream().mapToInt(downsampledCountFunGen.apply(rate)).sum();
                        } else {
                            return pool.submit(() ->
                                    umiCollections.parallelStream().mapToInt(downsampledCountFunGen.apply(rate))
                            ).get().sum();
                        }
                    } catch (InterruptedException | ExecutionException e) {
                        e.printStackTrace();
                    }
                    return null;
                })
                .collect(Collectors.toList());
    }

    /**
     * Given the number of transcripts found for a particular cell barcode and a map from quantiles (as an integer
     * between 0 and 100 representing percentage), returns the largest quantile q of
     * given that it has less than or equal to the median number of transcripts for that quantile.
     *
     * @param numTranscripts number of transcripts for a particular cell to be assigned to a quantile
     * @param quantiles map from different quantiles (Integer between 0 and 100 representing a percentage)
     *                  to the value of the median number of transcripts for cells in that quantile
     * @return the name of the quantile (Integer from 0 to 100) to which a given integer number of transcripts belongs
     */
    private Integer getQuantile(Integer numTranscripts, Map<Integer, Double> quantiles) {
        Integer[] keys = quantiles.keySet().stream().sorted().<Integer>toArray(Integer[]::new);
        for (int i = keys.length - 1; i >= 0; i--) {
            if (quantiles.get(keys[i]) <= numTranscripts)
                return keys[i];
        }
        return 0;
    }

    /**
     * Creates an object counter to keep track of how many cells are in each quantile.
     * Maps each quantile to a cumulative count (e.g. number of cells above the 30th percentile).
     *
     * @param transcriptsPerCell number of transcripts found for each distinct cell barcode
     * @param quantiles map of the number of transcripts per cell in each quantile
     * @return ObjectCounter of the number of cells that fall within the top i quantiles of numTranscripts per cell
     */
    private ObjectCounter<Integer> getCellQuantileCounts(ObjectCounter<String> transcriptsPerCell, Map<Integer, Double> quantiles) {
        ObjectCounter<Integer> res = new ObjectCounter<>();
        for (String cell : transcriptsPerCell.getKeys()) {
            Integer q = getQuantile(transcriptsPerCell.getCountForKey(cell), quantiles);
            quantiles.keySet().stream().filter(x->x<=q).forEach(x-> res.incrementByCount(x,1));
        }
        return res;
    }

    /**
     * Getter method for transcriptsPerCell
     * @return Object counter of transcripts per cell (equivalent to last column of OUTPUT_DOWNSAMPLING_FILE)
     */
    ObjectCounter<String> getTranscriptsPerCell() { return transcriptsPerCell; }

    /**
     * Method for running downsampling algorithm. Reads in INPUT and will print out one column per item (in order)
     * of DOWNSAMPLING_RATES to OUTPUT_DOWNSAMPLING_FILE.
     */
    private void writeDownsampledCellCounts() throws IOException {
        // Create output streams
        PrintStream writer = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(OUTPUT_DOWNSAMPLING_FILE));

        UMICollectionByCellParser iter = new UMICollectionByCellParser(INPUT);
        if (!iter.hasNext()) {
            throw new IOException("File "+INPUT+" is empty!");
        }
        List<UMICollection> cur;

        //Write header
        if (iter.hasNext()) {
            List<String> rates = DOWNSAMPLING_RATES.stream().map(String::valueOf).collect(Collectors.toList());
            Collections.replaceAll(rates, "1.0", "1");
            writer.println(StringUtil.join("\t", "CELL_BARCODE", String.join("\t", rates)));
        }

        int numThreads = (NUM_THREADS > 0? NUM_THREADS: Runtime.getRuntime().availableProcessors() + NUM_THREADS);
        pool = new ForkJoinPool(numThreads);
        while (iter.hasNext()) {
            // Get UMIs of next cell
            cur = iter.next();
            if (CELL_BC_FILE != null && !cellBarcodes.contains(cur.get(0).getCellBarcode()))
                continue;
            for (UMICollection c : cur) {
                transcriptsPerCell.incrementByCount(c.getCellBarcode(), c.getMolecularBarcodes().size());
            }
            writer.println(StringUtil.join("\t", cur.get(0).getCellBarcode(),
                    StringUtil.join("\t", getDownsampledCounts(cur))));
        }
        pool.shutdown();
        writer.flush();
        writer.close();
        CloserUtil.close(iter);
    }

    /**
     * Method for generating a summary of quantiles based on the number of transcripts for each cell.
     * To be called only after writeDownsampledCellCounts. Writes out to OUTPUT_QUANTILE_FILE
     */
    private void writeQuantileSummaryFile() {
        PrintStream writer = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(OUTPUT_QUANTILE_FILE));
        Map<Integer, Double> quantiles = percentiles().indexes(QUANTILES).compute(transcriptsPerCell.getCounts());
        ObjectCounter<Integer> quantileCounts = getCellQuantileCounts(transcriptsPerCell, quantiles);

        // Write header
        writer.println(StringUtil.join("\t", "quantile","cumulative_num_cells","median_transcripts"));

        // Largest cell num transcripts
        int largestCellSize = transcriptsPerCell.getCountForKey(transcriptsPerCell.getKeysOrderedByCount(true).get(0));
        writer.println(StringUtil.join("\t", "Largest_Cell", "1", String.valueOf(largestCellSize)));
        // Write out the rest except for 0; this should be the Smallest_Cell
        for (int i : Lists.reverse(QUANTILES)) {
            if (i == 100) continue;
            if (i == 0) break;
            writer.println(StringUtil.join("\t", (double)i/100, quantileCounts.getCountForKey(i), quantiles.get(i)));
        }
        writer.println(StringUtil.join("\t", "Smallest_Cell", quantileCounts.getCountForKey(0), quantiles.get(0)));
        writer.flush();
        writer.close();
    }

    @Override
    protected int doWork() {
        try {
            IOUtil.assertFileIsReadable(INPUT);
            IOUtil.assertFileIsWritable(OUTPUT_DOWNSAMPLING_FILE);
            IOUtil.assertFileIsWritable(OUTPUT_QUANTILE_FILE);
            if (CELL_BC_FILE != null)
                IOUtil.assertFileIsReadable(CELL_BC_FILE);
        } catch (Exception e) {
            e.printStackTrace();
            return 1;
        }
        if (CELL_BC_FILE != null)
            cellBarcodes = new HashSet<>(ParseBarcodeFile.readCellBarcodeFile(CELL_BC_FILE));
        if (RANDOM_SEED == null)
            random = new Random();
        else
            random = new Random(RANDOM_SEED);

        try {
            writeDownsampledCellCounts();
        } catch (IOException e) {
            e.printStackTrace();
            return 1;
        }
        writeQuantileSummaryFile();
        return 0;
    }

    public static void main(final String[] args) {
        new DownsampleTranscriptsAndQuantiles().instanceMainWithExit(args);
    }
}


