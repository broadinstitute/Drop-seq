/*
 * MIT License
 *
 * Copyright 2022 Broad Institute
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

import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.TabbedInputParser;
import picard.util.TabbedTextFileWithHeaderParser;
import picard.util.TabbedTextFileWithHeaderParser.Row;

import java.io.BufferedInputStream;
import java.io.File;
import java.util.HashMap;
import java.util.Map;

@CommandLineProgramProperties(
        summary =
                "Compare two sets of eQTLs counting the number of times the effect direction matches",
        oneLineSummary =
                "Compare two sets of eQTLs counting the number of times the effect direction matches",
        programGroup = DropSeq.class
)
public class SignTest extends CommandLineProgram {

    private final Log log = Log.getInstance(SignTest.class);

    @Argument(
            shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "The input eQTL file usually filtered to only contain the gene/variant pairs "
                    + "with the most significant effect."
    )
    public File INPUT;

    @Argument(
            shortName = "U",
            doc = "An unfiltered list of eQTLs with effects for a number of gene/variant combinations."
    )
    public File UNFILTERED_EQTL_FILE;

    @Argument(
            shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "The output sign test metrics."
    )
    public File OUTPUT;

    @Argument(
            shortName = "Q",
            doc = "Maximum q-value to include from the INPUT file. "
                    + "Set QVALUE_THRESHOLD=null to include all eQTLs.",
            optional = true
    )
    public Double QVALUE_THRESHOLD = 0.05;

    @Argument(
            shortName = "F",
            doc = "Also flip the reference and alternate alleles in the eQTLs from the UNFILTERED_EQTL_FILE that match the flipped "
                    + "alternate / reference alleles on the INPUT file variants."
    )
    public boolean FLIP_ALLELES = false;

    @Argument(
            shortName = "S",
            doc = "Strip Ensembl Gene suffix when comparing eQTLs"
    )
    public boolean STRIP_ENSG_SUFFIX = false;

    @Argument(doc = "Gene column label")
    public String GENE_COLUMN_LABEL = "phenotype_id";

    @Argument(doc = "Variant column label")
    public String VARIANT_COLUMN_LABEL = "variant_id";

    @Argument(doc = "Effect column label")
    public String EFFECT_COLUMN_LABEL = "slope";

    @Argument(doc = "q-value column label")
    public String QVALUE_COLUMN_LABEL = "qval";

    @Argument(doc = "Unfiltered gene column label")
    public String UNFILTERED_GENE_COLUMN_LABEL = "gene";

    @Argument(doc = "Unfiltered variant column label")
    public String UNFILTERED_VARIANT_COLUMN_LABEL = "SNP";

    @Argument(doc = "Unfiltered effect column label")
    public String UNFILTERED_EFFECT_COLUMN_LABEL = "beta";

    private static final int IO_BUFFER_SIZE = 10 * 1024 * 1024;

    public enum Effect {POSITIVE, NEGATIVE, UNKNOWN}

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(UNFILTERED_EQTL_FILE);
        IOUtil.assertFileIsWritable(OUTPUT);

        final SignTestMetric signTestMetric = new SignTestMetric();
        signTestMetric.FILE_UNFILTERED = UNFILTERED_EQTL_FILE.getAbsolutePath();
        signTestMetric.FILE_PERMUTED = INPUT.getAbsolutePath();
        signTestMetric.QVALUE_THRESHOLD = QVALUE_THRESHOLD;

        // For each gene, are the effect directions positive?
        final Map<String, Map<String, Effect>> inputGeneEffects = new HashMap<>();
        // For each gene, when a match is found, is the match's effect direction positive?
        final Map<String, Map<String, Effect>> matchGeneEffects = new HashMap<>();

        final TabbedTextFileWithHeaderParser inputParser =
                new TabbedTextFileWithHeaderParser(
                        new TabbedInputParser(
                                false,
                                new BufferedInputStream(
                                        IOUtil.openFileForReading(this.INPUT),
                                        IO_BUFFER_SIZE
                                )
                        )
                );

        log.info("Reading eQTLs from " + INPUT.getAbsolutePath());

        for (final Row inputRow : inputParser) {
            if (QVALUE_THRESHOLD != null) {
                try {
                    final String qvalueString = inputRow.getField(QVALUE_COLUMN_LABEL);
                    if (Double.parseDouble(qvalueString) > QVALUE_THRESHOLD) {
                        continue;
                    }
                } catch (final NumberFormatException e) {
                    continue;
                }
            }

            final String gene = stripEnsgSuffix(inputRow.getField(GENE_COLUMN_LABEL));
            final String variant = inputRow.getField(VARIANT_COLUMN_LABEL);
            final Effect effect = effectFromString(inputRow.getField(EFFECT_COLUMN_LABEL));

            inputGeneEffects.putIfAbsent(gene, new HashMap<>());
            inputGeneEffects.get(gene).putIfAbsent(variant, effect);

            matchGeneEffects.putIfAbsent(gene, new HashMap<>());

            signTestMetric.NUM_EQTLS_PERMUTED++;
        }

        inputParser.close();

        final TabbedTextFileWithHeaderParser unfilteredParser =
                new TabbedTextFileWithHeaderParser(
                        new TabbedInputParser(
                                false,
                                new BufferedInputStream(
                                        IOUtil.openFileForReading(this.UNFILTERED_EQTL_FILE),
                                        IO_BUFFER_SIZE
                                )
                        )
                );

        log.info("Reading unfiltered eQTLs from " + UNFILTERED_EQTL_FILE.getAbsolutePath());

        final ProgressLogger progress =
                new ProgressLogger(log, 10000000, "Processed", "unfiltered eQTLs");

        for (final Row unfilteredRow : unfilteredParser) {

            progress.record(null, 0);

            final String gene =
                    stripEnsgSuffix(unfilteredRow.getField(UNFILTERED_GENE_COLUMN_LABEL));
            final Map<String, Effect> inputVariantEffects = inputGeneEffects.get(gene);
            if (inputVariantEffects == null) {
                continue;
            }

            final Effect effect =
                    effectFromString(unfilteredRow.getField(UNFILTERED_EFFECT_COLUMN_LABEL));
            final Map<String, Effect> matchVariantEffects = matchGeneEffects.get(gene);

            final String variant = unfilteredRow.getField(UNFILTERED_VARIANT_COLUMN_LABEL);

            if (inputVariantEffects.containsKey(variant)) {
                matchVariantEffects.putIfAbsent(variant, effect);
            }

            /*
            Flipping here is slightly less efficient than flipping input eQTLs.
            However it's more consistent with the previous R implementation when processing
            unfiltered eQTLs that contain both the original and flipped indels.
            Ex: Both "C16orf87 / chr16:46811586:TA:T / +" AND "C16orf87 / chr16:46811586:T:TA / -"
             */
            if (FLIP_ALLELES) {
                final String flippedVariant = flipAlleles(variant);
                if (inputVariantEffects.containsKey(flippedVariant)) {
                    final Effect flippedEffect = flipEffect(effect);
                    matchVariantEffects.putIfAbsent(flippedVariant, flippedEffect);
                }
            }
        }

        unfilteredParser.close();

        log.info("Generating sign test results");

        for (final String gene : matchGeneEffects.keySet()) {
            final Map<String, Effect> matchVariantEffects = matchGeneEffects.get(gene);
            if (matchVariantEffects.isEmpty()) {
                continue;
            }
            final Map<String, Effect> inputVariantEffects = inputGeneEffects.get(gene);

            // For this gene are there matching gene/variant pairs, even ignoring effect direction?
            boolean geneHasMatchPair = false;
            // Do any of the matching pairs have matching effect directions?
            boolean geneAnyMatchEffect = false;
            // Do all of the matching pairs have matching effect directions?
            boolean geneAllMatchEffect = true;

            for (final String variant : matchVariantEffects.keySet()) {
                final Effect matchEffect = matchVariantEffects.get(variant);
                final Effect inputEffect = inputVariantEffects.get(variant);

                signTestMetric.NUM_EQTLS_COMPARED++;
                if (matchEffect != Effect.UNKNOWN) {
                    signTestMetric.NUM_EQTLS_UNFILTERED++;
                }

                geneHasMatchPair = true;
                if (matchEffect != Effect.UNKNOWN && matchEffect == inputEffect) {
                    geneAnyMatchEffect = true;
                    signTestMetric.NUM_EQTLS_SIGNS_MATCH++;
                } else {
                    geneAllMatchEffect = false;
                }
            }

            if (geneHasMatchPair) {
                signTestMetric.NUM_EGENES_UNIQUE++;
                if (geneAnyMatchEffect) {
                    signTestMetric.NUM_EGENES_ANY_MATCH++;
                    if (geneAllMatchEffect) {
                        signTestMetric.NUM_EGENES_ALL_MATCH++;
                    }
                }
            }
        }

        signTestMetric.FRAC_EQTLS_SIGNS_MATCH =
                signTestMetric.NUM_EQTLS_SIGNS_MATCH / (double) signTestMetric.NUM_EQTLS_COMPARED;
        signTestMetric.FRAC_EGENES_ANY_MATCH =
                signTestMetric.NUM_EGENES_ANY_MATCH / (double) signTestMetric.NUM_EGENES_UNIQUE;
        signTestMetric.FRAC_EGENES_ALL_MATCH =
                signTestMetric.NUM_EGENES_ALL_MATCH / (double) signTestMetric.NUM_EGENES_UNIQUE;

        log.info("Writing sign test results to " + OUTPUT.getAbsolutePath());

        final MetricsFile<SignTestMetric, ?> metricsFile = getMetricsFile();
        metricsFile.addMetric(signTestMetric);
        metricsFile.write(OUTPUT);

        return 0;
    }

    private String stripEnsgSuffix(final String gene) {
        if (STRIP_ENSG_SUFFIX && gene != null) {
            final int idx = gene.indexOf('.');
            return idx < 0 ? gene : gene.substring(0, idx);
        } else {
            return gene;
        }
    }

    private static Effect effectFromString(final String effect) {
        if (effect == null || effect.isEmpty()) {
            return Effect.UNKNOWN;
        }
        try {
            return Double.parseDouble(effect) > 0 ? Effect.POSITIVE : Effect.NEGATIVE;
        } catch (final NumberFormatException e) {
            return Effect.UNKNOWN;
        }
    }

    /**
     * Flip the effect direction for a slope / beta.
     */
    public static Effect flipEffect(final Effect effect) {
        //noinspection EnhancedSwitchMigration
        switch (effect) {
            case POSITIVE:
                return Effect.NEGATIVE;
            case NEGATIVE:
                return Effect.POSITIVE;
            default:
                return Effect.UNKNOWN;
        }
    }

    /**
     * Swap the ref and alt alleles.
     */
    private static String flipAlleles(final String variant) {
        final String[] variantTokens = variant.split(":", 4);
        final String chr = variantTokens[0];
        final String pos = variantTokens[1];
        final String ref = variantTokens[2];
        final String alt = variantTokens[3];
        return chr + ":" + pos + ":" + alt + ":" + ref;
    }

    public static class SignTestMetric extends MetricBase {
        /**
         * An unfiltered list of eQTLs with effects for a number of gene/variant combinations.
         */
        public String FILE_UNFILTERED;

        /**
         * The input eQTL file usually filtered to only contain the gene/variant pairs with the most significant effect.
         */
        public String FILE_PERMUTED;

        /**
         * Maximum q-value included from the FILE_PERMUTED file or null if all eQTLs were included.
         */
        public Double QVALUE_THRESHOLD;

        /**
         * Number of gene/variant pairs included from FILE_UNFILTERED.
         */
        public int NUM_EQTLS_UNFILTERED;

        /**
         * Number of gene/variant pairs included from FILE_PERMUTED.
         */
        public int NUM_EQTLS_PERMUTED;

        /**
         * How many gene/variant pairs have been located in both files?
         */
        public int NUM_EQTLS_COMPARED;

        /**
         * How many gene/variant pairs match effect directions?
         */
        public int NUM_EQTLS_SIGNS_MATCH;

        /**
         * NUM_EQTLS_SIGNS_MATCH / NUM_EQTLS_COMPARED
         */
        public double FRAC_EQTLS_SIGNS_MATCH;

        /**
         * How many genes have matching gene/variant pairs, even ignoring effect direction?
         */
        public int NUM_EGENES_UNIQUE;

        /**
         * How many genes have matching gene/variant pairs with a matching effect direction?
         */
        public int NUM_EGENES_ANY_MATCH;

        /**
         * How many genes have matching gene/variant pairs where all pairs match effect direction?
         */
        public int NUM_EGENES_ALL_MATCH;

        /**
         * NUM_EGENES_ANY_MATCH / NUM_EGENES_UNIQUE
         */
        public double FRAC_EGENES_ANY_MATCH;

        /**
         * NUM_EGENES_ALL_MATCH / NUM_EGENES_UNIQUE
         */
        public double FRAC_EGENES_ALL_MATCH;
    }

    /**
     * Stock main method.
     */
    public static void main(final String[] args) {
        System.exit(new SignTest().instanceMain(args));
    }
}
