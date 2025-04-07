/*
 * MIT License
 *
 * Copyright 2025 Broad Institute
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

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.SortingCollection;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import htsjdk.variant.vcf.VCFReader;
import htsjdk.variant.vcf.VCFRecordCodec;
import htsjdk.variant.vcf.VCFSimpleHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.OrderedConcurrentMapper;
import picard.PicardException;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.nio.PicardHtsPath;
import picard.util.TabbedInputParser;
import picard.util.TabbedTextFileWithHeaderParser;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.time.ZonedDateTime;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Objects;

@CommandLineProgramProperties(
        summary = "Convert a tab separated file of pair tests to a sites-only VCF file",
        oneLineSummary = "Convert a tab separated file of pair tests to a sites-only VCF file",
        programGroup = DropSeq.class
)
public class PairsToVcf extends CommandLineProgram {
    public PairsToVcf() {
        VALIDATION_STRINGENCY = ValidationStringency.LENIENT;
        CREATE_INDEX = true;
    }

    @Argument(
            shortName = StandardOptionDefinitions.INPUT_SHORT_NAME,
            doc = "Input TSV file."
    )
    public File INPUT;

    @Argument(
            shortName = "V",
            doc = "VCF file with the header and variants to use."
    )
    public PicardHtsPath VCF;

    @Argument(
            shortName = "M",
            doc = "A tab-separated file mapping old column names to new column names or VCF INFO fields.",
            optional = true
    )
    public File MAP_COLUMNS;

    @Argument(
            shortName = "C",
            doc = "Variant column in the TSV file."
    )
    public String VARIANT_COLUMN;

    @Argument(
            shortName = "S",
            doc = "Variant column field separator."
    )
    public String VARIANT_SEPARATOR = ":";

    @Argument(
            shortName = "B",
            doc = "Beta / Slope column in the TSV file.  If set the variants will be checked against the VCF " +
                    "and the beta value will be flipped if the variant genotype is swapped.",
            optional = true
    )
    public String BETA_COLUMN;

    @Argument(
            shortName = "T",
            doc = "Number of threads to use. Default is 1. If <1, the number of available processors is used.",
            optional = true
    )
    public Integer NUM_THREADS;

    @Argument(
            doc = "Set to false to disable sorting the output VCF file."
    )
    public Boolean SORT_OUTPUT = true;

    @Argument(
            doc = "Set to false to disable writing the command line to the VCF file."
    )
    public Boolean OUTPUT_COMMANDLINE = true;

    @Argument(
            shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output VCF file"
    )
    public File OUTPUT;

    private final Log log = Log.getInstance(PairsToVcf.class);

    private static final String SWAPPED_ANNOTATION = "VAS";
    private static final String SWAPPED_ANNOTATION_DESCRIPTION = "Variant alleles swapped";

    /**
     * How frequently to log progress.
     */
    private static final int PROGRESS_LOGGING_FREQUENCY = 1_000_000;

    /**
     * Items to cache so we don't have to re-query the VCF for the same variant multiple times.
     *
     * @param swapped If the variant information in the TSV is swapped relative to the VCF.
     * @param stop    The original stop if the variant genotype isn't swapped, otherwise the new stop position.
     */
    private record VariantContextInfo(boolean swapped, int stop) {
    }

    /**
     * Number of VariantContextInfo objects to cache. Useful since TSVs usually list the same variants multiple times.
     */
    private static final int VCI_CACHE_CAPACITY = 100_000;

    private final Map<String, VariantContextInfo> variantContextInfoCache = Collections.synchronizedMap(
            new LinkedHashMap<>(VCI_CACHE_CAPACITY, 0.75f, true) {
                @Override
                protected boolean removeEldestEntry(final Map.Entry<String, VariantContextInfo> eldest) {
                    return size() > VCI_CACHE_CAPACITY;
                }
            }
    );

    private Map<String, String> mapColumns;
    private List<String> columns;
    private SortingCollection<VariantContext> sorter;
    private VariantContextWriter writer;
    private ProgressLogger progressLogger;

    @Override
    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(VCF.toPath());
        IOUtil.assertFileIsWritable(OUTPUT);
        if (MAP_COLUMNS != null) {
            IOUtil.assertFileIsReadable(MAP_COLUMNS);
        }

        if (BETA_COLUMN != null && NUM_THREADS == null) {
            log.warn("BETA_COLUMN is set without NUM_THREADS, using default value of 1");
        }
        if (NUM_THREADS == null) {
            NUM_THREADS = 1;
        }

        readMapColumns();

        final SAMSequenceDictionary sequenceDictionary;
        try (final VCFFileReader reader = new VCFFileReader(VCF.toPath())) {
            sequenceDictionary = reader.getHeader().getSequenceDictionary();
        }

        try (final TabbedTextFileWithHeaderParser tsvParser = newTsvParser();
             final VariantContextWriter writer = newOutputWriter(sequenceDictionary)) {

            this.columns = tsvParser.columnLabelsList();

            final VCFHeader outputHeader = newVcfHeader(sequenceDictionary);
            writer.writeHeader(outputHeader);
            this.writer = writer;

            if (SORT_OUTPUT) {
                this.sorter = newOutputSorter(outputHeader);
            }
            final String loopVerb = SORT_OUTPUT ? "converted" : "wrote";
            this.progressLogger = new ProgressLogger(log, PROGRESS_LOGGING_FREQUENCY, loopVerb, "pairs");

            new OrderedConcurrentMapper
                    .Builder<TabbedTextFileWithHeaderParser.Row, VCFReader, VariantContext>()
                    .withNumThreads(NUM_THREADS)
                    .withMaxItems(MAX_RECORDS_IN_RAM)
                    .withIterable(tsvParser)
                    .withPerThreadFactory(this::newVcfReader)
                    .withMapper(this::mapRow)
                    .withConsumer(this::addOutputVariant)
                    .build()
                    .execute();

            writeSortedOutput();

        } finally {
            if (this.sorter != null) {
                this.sorter.cleanup();
            }
        }
        return 0;
    }

    /**
     * Read the MAP_COLUMNS file if it exists.
     */
    private void readMapColumns() {
        if (MAP_COLUMNS == null) {
            this.mapColumns = Collections.emptyMap();
            return;
        }

        final Map<String, String> mapColumns = new HashMap<>();
        try (final BufferedReader reader = new BufferedReader(new FileReader(MAP_COLUMNS))) {
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.strip();
                if (line.isEmpty() || line.startsWith("#")) {
                    continue;
                }
                final String[] parts = line.split("\t");
                mapColumns.put(parts[0], parts[1]);
            }
        } catch (IOException e) {
            throw new PicardException("Error processing MAP_COLUMNS file " + MAP_COLUMNS, e);
        }
        this.mapColumns = Collections.unmodifiableMap(mapColumns);
    }

    /**
     * Returns a buffered TSV parser.
     *
     * @return the TSV parser
     */
    private TabbedTextFileWithHeaderParser newTsvParser() {
        final InputStream fis = IOUtil.openFileForReading(INPUT);
        final BufferedInputStream inputStream = new BufferedInputStream(fis, 1024 * 1024);
        return new TabbedTextFileWithHeaderParser(new TabbedInputParser(false, inputStream));
    }

    /**
     * Create a new asynchronous VCF writer without indexing on the fly and without writing genotypes.
     *
     * @param sequenceDictionary The contigs that should match the TSV file
     * @return The VCF writer
     */
    private VariantContextWriter newOutputWriter(final SAMSequenceDictionary sequenceDictionary) {
        final VariantContextWriterBuilder builder = new VariantContextWriterBuilder();
        builder.setOutputFile(OUTPUT);
        builder.setReferenceDictionary(sequenceDictionary);
        builder.setOption(Options.DO_NOT_WRITE_GENOTYPES);
        builder.setOption(Options.USE_ASYNC_IO);
        builder.setBuffer(8192);
        builder.setCreateMD5(this.CREATE_MD5_FILE);
        if (this.CREATE_INDEX) {
            builder.setOption(Options.INDEX_ON_THE_FLY);
        } else {
            builder.unsetOption(Options.INDEX_ON_THE_FLY);
        }
        return builder.build();
    }

    /**
     * Create a new sorter for the output VCF file.
     *
     * @param outputHeader The VCF header
     * @return The sorting collection
     */
    private SortingCollection<VariantContext> newOutputSorter(final VCFHeader outputHeader) {
        final VCFRecordCodec codec = new VCFRecordCodec(
                outputHeader,
                this.VALIDATION_STRINGENCY != ValidationStringency.STRICT
        );
        return SortingCollection.newInstanceFromPaths(
                VariantContext.class,
                codec,
                outputHeader.getVCFRecordComparator(),
                this.MAX_RECORDS_IN_RAM,
                this.TMP_DIR.stream().map(File::toPath).toList()
        );
    }

    /**
     * Generate the VCF header including the tool name, version, date, command line, and column headers.
     *
     * @param sequenceDictionary The contigs to write to the VCF header
     * @return The VCF header
     */
    private VCFHeader newVcfHeader(final SAMSequenceDictionary sequenceDictionary) {
        final VCFHeader vcfHeader = new VCFHeader();
        final String toolName = this.getClass().getSimpleName();
        vcfHeader.addMetaDataLine(new VCFHeaderLine("source", toolName));

        if (OUTPUT_COMMANDLINE) {
            final Map<String, String> simpleHeaderLineMap = new HashMap<>(4);
            String version = this.getVersion();
            if (version == null) {
                version = "unknown";
            } else if (version.startsWith("Version:")) {
                version = version.substring("Version:".length());
            }
            simpleHeaderLineMap.put("ID", toolName);
            simpleHeaderLineMap.put("Version", version);
            simpleHeaderLineMap.put("Date", ZonedDateTime.now().toString());
            simpleHeaderLineMap.put("CommandLine", this.getCommandLine());
            vcfHeader.addMetaDataLine(new VCFSimpleHeaderLine(toolName + "CommandLine", simpleHeaderLineMap));
        }

        vcfHeader.setSequenceDictionary(sequenceDictionary);

        if (BETA_COLUMN != null) {
            vcfHeader.addMetaDataLine(new VCFFilterHeaderLine(SWAPPED_ANNOTATION, SWAPPED_ANNOTATION_DESCRIPTION));
        }

        this.columns.forEach(column ->
                {
                    if (Objects.equals(column, VARIANT_COLUMN)) return;
                    final String mappedColumn = mapColumns.getOrDefault(column, column);
                    final VCFInfoHeaderLine headerLine =
                            new VCFInfoHeaderLine(mappedColumn, 1, VCFHeaderLineType.String, column);
                    vcfHeader.addMetaDataLine(headerLine);
                }
        );

        return vcfHeader;
    }

    private VCFReader newVcfReader() {
        return BETA_COLUMN == null ? null : new VCFFileReader(VCF.toPath());
    }

    /**
     * Map a TSV row to a VCF record.
     *
     * @param tsvRow    The variant row
     * @param vcfReader The thread local VCF reader instance
     * @return The VCF record with the TSV row information
     */
    private VariantContext mapRow(final TabbedTextFileWithHeaderParser.Row tsvRow, final VCFReader vcfReader) {
        final String variant = tsvRow.getField(VARIANT_COLUMN);
        final String[] parts = variant.split(VARIANT_SEPARATOR);
        final String contig = parts[0];
        final int position = Integer.parseInt(parts[1]);
        final String refVariant = parts[2];
        final String altVariant = parts[3];

        final VariantContextBuilder builder = new VariantContextBuilder();
        builder.chr(contig);
        builder.start(position);
        if (BETA_COLUMN == null) {
            builder.stop(position + refVariant.length() - 1);
            builder.alleles(refVariant, altVariant);
        } else {
            final VariantContextInfo result = query(vcfReader, variant);
            builder.stop(result.stop());

            final String betaColumn = mapColumns.getOrDefault(BETA_COLUMN, BETA_COLUMN);
            final String betaStr = tsvRow.getField(BETA_COLUMN);
            if (result.swapped()) {
                builder.filter(SWAPPED_ANNOTATION);
                builder.alleles(altVariant, refVariant);
                try {
                    final double betaDbl = Double.parseDouble(betaStr);
                    // avoid -0
                    builder.attribute(betaColumn, betaDbl < 0 || betaDbl > 0 ? String.valueOf(-betaDbl) : "0");
                } catch (NumberFormatException e) {
                    // Anything non-numerical is left as is
                    builder.attribute(betaColumn, betaStr);
                }
            } else {
                builder.alleles(refVariant, altVariant);
                builder.attribute(betaColumn, betaStr);
            }
        }

        for (final String column : this.columns) {
            if (!column.equals(VARIANT_COLUMN) && !column.equals(BETA_COLUMN)) {
                final String mappedColumn = mapColumns.getOrDefault(column, column);
                final String value = tsvRow.getField(column);
                builder.attribute(mappedColumn, value);
            }
        }
        return builder.make();
    }

    /**
     * Find the VCF record for the given variant.
     *
     * <p>This operation is relatively expensive as it needs to both seek and decompress the VCF records.
     * So we cache the results.
     *
     * @param vcfReader Thread local VCF reader instance
     * @param variant   The variant to find
     * @return The variant context information
     */
    private VariantContextInfo query(final VCFReader vcfReader, final String variant) {
        VariantContextInfo result = variantContextInfoCache.get(variant);
        if (result != null) {
            return result;
        }

        final String[] parts = variant.split(VARIANT_SEPARATOR);
        final String contig = parts[0];
        final int position = Integer.parseInt(parts[1]);
        final String refVariant = parts[2];
        final String altVariant = parts[3];
        try (final CloseableIterator<VariantContext> query = vcfReader.query(contig, position, position)) {
            while (query.hasNext()) {
                final VariantContext next = query.next();
                /*
                Skip if the position does not match, for example two overlapping indels at different positions.
                #CHROM  POS       ID  REF     ALT
                chr4    41086938  .   ACACAC  A
                chr4    41086940  .   ACACAC  A
                 */
                if (position != next.getStart()) {
                    continue;
                }
                final String refVcf = next.getReference().getBaseString();
                for (final Allele altAllele : next.getAlternateAlleles()) {
                    final String altVcf = altAllele.getBaseString();
                    if (refVariant.equals(refVcf) && altVariant.equals(altVcf)) {
                        result = new VariantContextInfo(false, next.getEnd());
                        break;
                    } else if (result == null && refVariant.equals(altVcf) && altVariant.equals(refVcf)) {
                        result = new VariantContextInfo(true, position + altVariant.length() - 1);
                    }
                }
                if (result != null && !result.swapped()) {
                    break;
                }
            }
        }
        if (result == null) {
            switch (VALIDATION_STRINGENCY) {
                case STRICT:
                    throw new IllegalArgumentException("VCF does not contain " + variant);
                case LENIENT:
                    log.warn("VCF does not contain " + variant); // then fall through
                case SILENT:
                    result = new VariantContextInfo(false, position + refVariant.length() - 1);
            }
        }
        variantContextInfoCache.put(variant, result);
        return result;
    }

    /**
     * Output a variant created from a TSV row.
     *
     * @param vc The variant to output
     */
    private void addOutputVariant(final VariantContext vc) {
        if (this.sorter != null) {
            this.sorter.add(vc);
        } else {
            this.writer.add(vc);
        }
        this.progressLogger.record(vc.getContig(), vc.getStart());
    }

    private void writeSortedOutput() {
        if (this.sorter != null) {
            final ProgressLogger writingLogger =
                    new ProgressLogger(log, PROGRESS_LOGGING_FREQUENCY, "wrote", "pairs");
            try (final CloseableIterator<VariantContext> iterator = this.sorter.iterator()) {
                while (iterator.hasNext()) {
                    final VariantContext next = iterator.next();
                    this.writer.add(next);
                    writingLogger.record(next.getContig(), next.getStart());
                }
            }
        }
    }

    public static void main(final String[] args) {
        new PairsToVcf().instanceMainWithExit(args);
    }
}
