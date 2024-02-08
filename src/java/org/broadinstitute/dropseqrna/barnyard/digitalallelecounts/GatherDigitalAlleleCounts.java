/*
 * MIT License
 *
 * Copyright 2017 Broad Institute
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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import java.io.File;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.collections4.SetUtils;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.BarcodeListRetrieval;
import org.broadinstitute.dropseqrna.barnyard.GeneFunctionCommandLineBase;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.CountChangingIteratorWrapper;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.VCFUtils;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.GeneFunctionIteratorWrapper;
import org.broadinstitute.dropseqrna.utils.readiterators.IgnoreGeneAnnotationTagger;
import org.broadinstitute.dropseqrna.utils.readiterators.PCRDuplicateFilteringIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.statistics.BinomialStatistics;
import org.broadinstitute.dropseqrna.vcftools.SampleAssignmentVCFUtils;
import org.broadinstitute.dropseqrna.vcftools.filters.HetSNPFilter;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import picard.annotation.LocusFunction;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(summary = "Measures the digital allele counts of a library.  " + "Method: 1) Order BAM by SNP/GENE/Cell/UMI tags."
		+ "2) Build per snp/gene/cell/umi pileups of reads" + "3) Filter pileups by base quality, count reads on each UMI"
		+ "4) Gather digital allele counts on a SNP/GENE/CELL" + "5) Repeat across all cells of the SNP/GENE" + "6) Repeat across all SNPs/GENEs"
		+ "This program requires a tag for what gene a read is on, a molecular barcode tag, and a exon tag.  The exon and gene tags may not be present on every read."
		+ "When filtering the data for a set of barcodes to use, the data is filtered by ONE of the following methods (and if multiple params are filled in, the top one takes precidence):"
		+ "1) CELL_BC_FILE, to filter by the some fixed list of cell barcodes" + "2) MIN_NUM_GENES_PER_CELL " + "3) MIN_NUM_TRANSCRIPTS_PER_CELL "
		+ "4) NUM_CORE_BARCODES " + "5) MIN_NUM_READS_PER_CELL", oneLineSummary = "Calculate Digital Allele Counts", programGroup = DropSeq.class)

public class GatherDigitalAlleleCounts extends GeneFunctionCommandLineBase {
	private static final Log log = Log.getInstance(GatherDigitalAlleleCounts.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze. This argument can accept wildcards, or a file with the suffix .bam_list that contains the locations of multiple BAM files", minElements = 1)
	public List<File> INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output file of DAC Results.  One SNP/GENE/CELL per row. This supports zipped formats like gz and bz2.", optional = true)
	public File OUTPUT;

	@Argument(doc = "A file with 2 columns, defining how cells should be grouped for meta analysis.  The first column contains the cluster identifier, the second the cell barcode.  This file is tab separated, and has a header with the column names CLUSTER BARCODE", optional = true)
	public File CLUSTER_FILE;

	@Argument(doc = "If a cluster file is present to cluster cells into meta-analysis results, this parameter directs the output to the given file.", optional = true)
	public File CLUSTER_OUTPUT;

	@Argument(doc = "Output file of DAC meta analysis results.  One SNP/GENE per row. The results of all cells are combined, and joint pvalues and confidence thresholds calculated.  This supports zipped formats like gz and bz2.", optional = true)
	public File SUMMARY;

	@Argument(doc = "Output file containing the SNP pileup aggregated across all cells/genes, and the inferred allele frequency.", optional = true)
	public File ALLELE_FREQUENCY_OUTPUT;

	/**
	 * VCF OPTIONS
	 */

	@Argument(doc = "The input VCF file to analyze.")
	public File VCF;

	@Argument(doc = "A file with a list of samples in the VCF.  This subsets the VCF into a smaller data set containing only the samples listed. The file has 1 column with no header.", optional = true)
	public File SAMPLE_FILE = null;

	@Argument(doc = "Directly set the list of samples to reference in the VCF.  This option can be specified multiple times.  This option superceeds SAMPLE_FILE.", optional = true)
	public List<String> SAMPLEs;

	@Argument(doc = "The minimum genotype quality for a genotype to be counted as called for a donor.")
	public Integer GQ_THRESHOLD = 30;

	@Argument(doc = "At least <FRACTION_SAMPLES_PASSING> samples must have genotype scores >= GQ_THRESHOLD for the variant in the VCF to be included in the analysis.")
	public double FRACTION_SAMPLES_PASSING = 0.5;

	@Argument(doc = "A list of chromosomes to omit from the analysis.  The default is to omit the sex chromosomes.")
	public List<String> IGNORED_CHROMOSOMES = null;

	@Argument(doc = "Only output heterozygous SNPs for these samples.  All samples must be heterozygous for the SNP to be analyzed.")
	public boolean HET_SNPS_ONLY = false;

	@Argument(doc = "Only use SNPs that vary in the population of samples to analyze in the VCF")
	public boolean POLYMORPHIC_SNPS_ONLY = true;

	/***
	 * Processing options
	 */

	@Argument(doc = "In cases where there are multiple variants per read, the read can be counted multiple times.  If set to true, only counts one variant per read."
			+ "The variant selected is the one with the highest mean genotype quality.")
	public boolean SINGLE_VARIANT_READS = false;

	@Argument(doc = "In cases where a read overlaps a gene, set this to true to count that SNP on both genes.  When set to false, the genes are sorted by the following criteria:"
			+ "UMI purity, number of UMIs, average base quality of the reads, and number of supporting reads.  The best gene is selected from a list of possible genes."
			+ "This should be set to false when using the ALLELE_FREQUENCY_OUTPUT to avoid repeated SNP positions.")
	public boolean MULTI_GENES_PER_READ = true;

	@Argument(doc = "[EXPERIMENTAL] Special case adaptation for reads that are not tagged with gene names and functions.  These reads are typically INTERGENIC reads [which don't have genes by definition "
			+ "but can also come from protocols that have UMIs but do not make sense to align to genes.  Gene names for these reads (if any) will be replaced by the contig name.  "
			+ "Tagged strand will be oded to the read strand, and gene function will be set to first locus function in the LOCUS_FUNCTION_LIST.  Pileups will be generated by all "
			+ "reads that intersect each SNP pileup. This has the added benefit of accessing intergenic regions, which normally do not have genes.")
	public boolean INCLUDE_UNTAGGED_READS = false;

	@Argument(doc = "[EXPERIMENTAL] If INCLUDE_UNTAGGED_READS is true, then this option splits the contig results into the positive and negative strand.")
	public boolean SPLIT_CONTIG_BY_STRAND = false;

	/***
	 * Tags and filtering options.
	 */

	@Argument(doc = "The cell barcode tag.  If there are no reads with this tag, the program will assume that all reads belong to the same cell and process in single sample mode.")
	public String CELL_BARCODE_TAG = "XC";

	@Argument(doc = "The molecular barcode tag.")
	public String MOLECULAR_BARCODE_TAG = "XM";

	@Argument(doc = "The edit distance that molecular barcodes should be combined at within a gene/SNP.")
	public Integer EDIT_DISTANCE = 1;

	@Argument(doc = "The map quality of the read to be included.")
	public Integer READ_MQ = 10;

	@Argument(doc = "The base quality of the read bases to be included.")
	public Integer BASE_QUALITY = 10;

	@Argument(doc = "The confidence score used to calculate the binomial confidence intervals.")
	public Double CONFIDENCE = 0.95;

	@Argument(doc = "Gather up all cell barcodes that have more than some number of reads.", optional = true)
	public Integer MIN_NUM_READS_PER_CELL = null;

	@Argument(doc = "Number of cells that you think are in the library.  This accomplishes the same goals as the MIN_NUM_READS_CORE argument, but instead of defining barcodes as important based on the number of reads, it picks the top <X> barcodes as core.", optional = true)
	public Integer NUM_CORE_BARCODES = null;

	@Argument(doc = "Override CELL_BARCODE and MIN_NUM_READS_PER_CELL, and process reads that have the cell barcodes in this file instead.  The file has 1 column with no header.", optional = true)
	public File CELL_BC_FILE = null;

	@Argument(doc = "Remove UMIs below this purity threshold.  A UMI's purity is determined as the number of reads of the most common base divided by the total number of reads.")
	public Double UMI_PURITY_THRESHOLD = 1.0;

	@Argument(doc = "Skip duplicate reads.")
	public boolean FILTER_PCR_DUPES = true;

	// private final boolean AUTO_FLUSH_OUTPUTS=true;
	private final String SNP_TAG = "YS";
	private final char[] BASES = { 'A', 'C', 'G', 'T', 'N' };
	private final int PROGRESS_RATE = 100000;

	@Override
	protected String[] customCommandLineValidation() {

		final ArrayList<String> list = new ArrayList<>(1);
		IOUtil.assertFileIsReadable(this.VCF);
		if (this.SAMPLE_FILE != null)
			IOUtil.assertFileIsReadable(this.SAMPLE_FILE);
		if (this.CLUSTER_FILE != null) {
			if (this.CLUSTER_OUTPUT == null)
				throw new IllegalArgumentException("If the CLUSTER_FILE parameter is set, the CLUSTER_OUTPUT must also be set.");
			IOUtil.assertFileIsReadable(this.CLUSTER_FILE);
			IOUtil.assertFileIsWritable(this.CLUSTER_OUTPUT);
		}
		if (this.OUTPUT != null)
			IOUtil.assertFileIsWritable(this.OUTPUT);
		if (this.SUMMARY != null)
			IOUtil.assertFileIsWritable(this.SUMMARY);
		if (this.ALLELE_FREQUENCY_OUTPUT != null)
			IOUtil.assertFileIsWritable(this.ALLELE_FREQUENCY_OUTPUT);

		if (INCLUDE_UNTAGGED_READS == false && SPLIT_CONTIG_BY_STRAND)
			list.add("If INCLUDE_UNTAGGED_READS is false, can not set SPLIT_CONTIG_BY_STRAND to true!");

		// TODO: at some point couple this more tightly to IgnoreGeneAnnotationTagger. This message will be wrong if that
		// implementation changes.
		if (INCLUDE_UNTAGGED_READS)
			log.info("INCLUDE_UNTAGGED_READS is true.  Setting untagged reads to have gene function [" + this.LOCUS_FUNCTION_LIST.get(0).toString());

		// VALIDATE that either CELL_BC_FILE or NUM_CORE_BARCODES is not null
		if (this.CELL_BC_FILE == null && this.NUM_CORE_BARCODES == null)
			list.add("Must set either CELL_BC_FILE or NUM_CORE_BARCODES argument!");

		return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);

	}

	@Override
	protected int doWork() {
		this.INPUT = FileListParsingUtils.expandFileList(INPUT);

		PrintStream outCluster = null;
		Map<String, Set<String>> clusterMap = null;
		// if CLUSTER_FILE is set and the output is set, prepare the outputs.
		if (this.CLUSTER_FILE != null && this.CLUSTER_OUTPUT != null) {
			outCluster = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(CLUSTER_OUTPUT));
			clusterMap = ParseBarcodeFile.readCellClusterFile(this.CLUSTER_FILE);
			writeClusterHeader(outCluster);
		}

		PrintStream out = getPrintStreamOrNull(this.OUTPUT);
		if (out != null)
			writeHeader(out);

		PrintStream outSum = getPrintStreamOrNull(this.SUMMARY);
		if (outSum != null)
			writeHeader(outSum);

		GdacAlleleFrequencyWriter afWriter = null;
		if (this.ALLELE_FREQUENCY_OUTPUT != null) {
			afWriter = new GdacAlleleFrequencyWriter(this.ALLELE_FREQUENCY_OUTPUT);
			afWriter.writeHeader();
		}

		List<String> cellBarcodes = new BarcodeListRetrieval().getCellBarcodes(this.INPUT, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG,
				this.GENE_NAME_TAG, this.GENE_STRAND_TAG, this.GENE_FUNCTION_TAG, this.STRAND_STRATEGY, this.LOCUS_FUNCTION_LIST, this.FUNCTIONAL_STRATEGY, this.CELL_BC_FILE,
				this.READ_MQ, null, null, this.MIN_NUM_READS_PER_CELL, this.NUM_CORE_BARCODES, this.EDIT_DISTANCE, null);

		log.info("Selected cell barcodes for analysis [" + cellBarcodes.size() + "]");

		SamHeaderAndIterator headerAndIter = SamFileMergeUtil.mergeInputs(this.INPUT, false, SamReaderFactory.makeDefault());

		SNPInfoCollection snpInfo = getSNPInfoCollection();

		if (snpInfo.getIntervalList().size() == 0) {
			log.error("No SNPs found, exiting!");
			return 1;
		}

		Map<Interval, Double> genotypeQuality = null;
		// select only the highest GQ SNP per read when enabled.
		if (this.SINGLE_VARIANT_READS)
			genotypeQuality = snpInfo.getAverageGQ();

		// override the normal gene annotations with new ones before any other operations.
		if (this.INCLUDE_UNTAGGED_READS) {
			// This only overwrites tags where the read does not match an accepted locus function.
			IgnoreGeneAnnotationTagger tagger = new IgnoreGeneAnnotationTagger(headerAndIter.iterator, this.GENE_NAME_TAG, this.GENE_STRAND_TAG,
					this.GENE_FUNCTION_TAG, this.LOCUS_FUNCTION_LIST, this.SPLIT_CONTIG_BY_STRAND, false);
			headerAndIter = new SamHeaderAndIterator(headerAndIter.header, (CloseableIterator<SAMRecord>) tagger.iterator());
		}

		if (this.FILTER_PCR_DUPES) {
			final PCRDuplicateFilteringIterator pcrDuplicateFilteringIterator = new PCRDuplicateFilteringIterator(headerAndIter.iterator);
			headerAndIter = new SamHeaderAndIterator(headerAndIter.header, pcrDuplicateFilteringIterator);
		}

		//TODO: should this be assignReadsToAllGenes be set to false to mirror DigitalExpression and donor assignment code?
		SNPUMIBasePileupIterator sbpi = new SNPUMIBasePileupIterator(headerAndIter, snpInfo.getIntervalList(), GENE_NAME_TAG, GENE_STRAND_TAG,
				GENE_FUNCTION_TAG, LOCUS_FUNCTION_LIST, STRAND_STRATEGY, this.FUNCTIONAL_STRATEGY, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, this.SNP_TAG, null, this.READ_MQ,
				true, cellBarcodes, genotypeQuality, SortOrder.SNP_GENE);

		MultiCellDigitalAlleleCountsIterator multiIter = getDACIterator(sbpi, MULTI_GENES_PER_READ, snpInfo);

		// sort cell barcodes alphabetically for output.
		Collections.sort(cellBarcodes);
		@SuppressWarnings("unused")
		int counter = 0;
		while (multiIter.hasNext()) {
			MultiCellDigitalAlleleCounts mcdac = multiIter.next();
			processMCDAC(cellBarcodes, mcdac, out, outSum, outCluster, afWriter, clusterMap);

			counter++;
			if (counter % PROGRESS_RATE == 0)
				log.info("Processed " + counter + " SNPs");
		}
		log.info("Processed " + counter + " total SNPs");
		CloserUtil.close(out);
		CloserUtil.close(outSum);
		CloserUtil.close(outCluster);
		CloserUtil.close(afWriter);
		multiIter.close();
		return 0;
	}

	/**
	 * Takes all the cells for a SNP/GENE and writes out their summary statistics.
	 * 
	 * @param cellBarcodes The list of cells to try and write out - this provides output ordering once you're past the
	 *                     SNP/GENE level.
	 * @param mcdac
	 * @param out          The per-cell output stream
	 * @param outSum       The meta analysis output stream
	 */
	private void processMCDAC(final List<String> cellBarcodes, final MultiCellDigitalAlleleCounts mcdac, final PrintStream out, final PrintStream outSum,
			final PrintStream outCluster, final GdacAlleleFrequencyWriter afOut, final Map<String, Set<String>> clusterMap) {

		// collapses all DAC UMIs.
		mcdac.collapseDACs(this.EDIT_DISTANCE);
		mcdac.filterDataByUMIPurity(UMI_PURITY_THRESHOLD);
		// write meta analysis stats for non-empty objects.
		DigitalAlleleCounts metaDAC = mcdac.getMetaAnalysis();
		// if the meta analysis has no information (all reads filtered out?) return.
		if (metaDAC.isEmpty())
			return;

		if (!metaDAC.isEmpty()) {
			writeStats(metaDAC, this.CONFIDENCE, outSum);
			if (afOut != null)
				afOut.writeLine(metaDAC);
			// if (autoFlush & outSum!=null) outSum.flush();
		}

		// if your cluster map and output writer are not null.
		if (outCluster != null && clusterMap != null) {
			// loop over clusters and write stats for non-empty objects.
			List<String> clusterIDs = new ArrayList<>(clusterMap.keySet());
			Collections.sort(clusterIDs);

			for (String clusterID : clusterIDs) {
				Set<String> barcodes = clusterMap.get(clusterID);
				DigitalAlleleCounts dac = mcdac.getMetaAnalysis(clusterID, barcodes);
				if (dac != null && !dac.isEmpty()) {
					writeStats(dac, this.CONFIDENCE, outCluster);
					// if (autoFlush) outCluster.flush();
				}
			}
		}

		// loop over cells and write stats for non-empty objects.
		for (String cell : cellBarcodes) {
			DigitalAlleleCounts dac = mcdac.getDigitalAlleleCounts(cell);
			// log.info(dac.getSnpInterval() + " " +dac.getCell() + " " +dac.getGene());

			if (dac != null && !dac.isEmpty()) {
				writeStats(dac, this.CONFIDENCE, out);
				// if (autoFlush) out.flush();
			}
		}
	}

	private void writeStats(final DigitalAlleleCounts dac, final double confidence, final PrintStream out) {
		if (out == null)
			return;
		DecimalFormat ratioFormat = new DecimalFormat("0.000");
		DecimalFormat pvalueFormat = new DecimalFormat("0.####E0");

		List<String> line = new ArrayList<>();
		// chromosome, position, gene, cell, reads, umis, numUMIs, read ratio, read pval, read CI, umi raio, umi pval, umi CI
		line.add(dac.getSnpInterval().getContig());
		line.add(Integer.toString(dac.getSnpInterval().getStart()));
		line.add(dac.getGene());
		line.add(dac.getCell());
		line.add(dac.getReferenceAllele() + "");
		line.add(dac.getAltAllele() + "");

		ObjectCounter<Character> readCounts = dac.getReadCounts();
		for (Character base : BASES)
			line.add(Integer.toString(readCounts.getCountForKey(base)));

		ObjectCounter<Character> umiCounts = dac.getUMIAlleleCount();
		for (Character base : BASES)
			line.add(Integer.toString(umiCounts.getCountForKey(base)));

		line.add(Integer.toString(dac.umis().size()));
		line.add(ratioFormat.format(dac.getMeanUMIPurity()));

		// if the number of trials is 0, this causes issues.
		int r = readCounts.getCountForKey(dac.getReferenceAllele()) + readCounts.getCountForKey(dac.getAltAllele());
		if (r > 0) {
			BinomialStatistics stats = dac.getBinomialStatistics(readCounts, confidence);
			line.add(pvalueFormat.format(stats.getBinomialPvalue()));
			line.add(ratioFormat.format(stats.getRatio()));
			line.add(ratioFormat.format(stats.getBinomialConfidenceInterval().getLowerBound()));
			line.add(ratioFormat.format(stats.getBinomialConfidenceInterval().getUpperBound()));
		} else {
			line.addAll(Arrays.asList("NA", "NA", "NA", "NA"));
		}

		r = umiCounts.getCountForKey(dac.getReferenceAllele()) + readCounts.getCountForKey(dac.getAltAllele());
		if (r > 0) {
			BinomialStatistics stats = dac.getBinomialStatistics(umiCounts, confidence);
			line.add(pvalueFormat.format(stats.getBinomialPvalue()));
			line.add(ratioFormat.format(stats.getRatio()));
			line.add(ratioFormat.format(stats.getBinomialConfidenceInterval().getLowerBound()));
			line.add(ratioFormat.format(stats.getBinomialConfidenceInterval().getUpperBound()));
		} else {
			line.addAll(Arrays.asList("NA", "NA", "NA", "NA"));
		}

		String h = StringUtils.join(line, "\t");
		out.println(h);
	}

	private void writeHeader(final PrintStream out) {
		if (out == null)
			return;
		// chromosome, position, gene, cell, reads, umis, numUMIs, read ratio, read pval, read CI, umi raio, umi pval, umi CI
		List<String> header = new ArrayList<>();

		header.add("chr");
		header.add("pos");
		header.add("gene");
		header.add("cell");
		header.add("ref_base");
		header.add("alt_base");
		for (char c : BASES)
			header.add("r" + c);
		for (char c : BASES)
			header.add(c + "");

		header.add("num_umi");
		header.add("umi_mean_purity");

		header.add("read_pval");
		header.add("read_ratio");
		header.add("read_ci_low");
		header.add("read_ci_high");

		header.add("umi_pval");
		header.add("umi_ratio");
		header.add("umi_ci_low");
		header.add("umi_ci_high");

		String h = StringUtils.join(header, "\t");
		out.println(h);
	}

	private void writeClusterHeader(final PrintStream out) {
		if (out == null)
			return;
		// chromosome, position, gene, cell, reads, umis, numUMIs, read ratio, read pval, read CI, umi raio, umi pval, umi CI
		List<String> header = new ArrayList<>();

		header.add("chr");
		header.add("pos");
		header.add("gene");
		header.add("clusterID");
		header.add("ref_base");
		for (char c : BASES)
			header.add("r" + c);
		for (char c : BASES)
			header.add(c + "");

		header.add("num_umi");
		header.add("umi_mean_purity");

		header.add("read_pval");
		header.add("read_ratio");
		header.add("read_ci_low");
		header.add("read_ci_high");

		header.add("umi_pval");
		header.add("umi_ratio");
		header.add("umi_ci_low");
		header.add("umi_ci_high");

		String h = StringUtils.join(header, "\t");
		out.println(h);
	}

	private SNPInfoCollection getSNPInfoCollection() {
		final VCFFileReader vcfReader = new VCFFileReader(this.VCF, false);

		if (!VCFUtils.GQInHeader(vcfReader)) {
			this.GQ_THRESHOLD = -1;
			log.info("Genotype Quality [GQ] not found in header.  Disabling GQ_THRESHOLD parameter");
		}

		List<String> vcfSamples = this.SAMPLEs;
		// if the vcfSamples isn't set, use the sample filee
		if (vcfSamples == null || vcfSamples.size() == 0)
			vcfSamples = SampleAssignmentVCFUtils.getVCFSamples(vcfReader, this.SAMPLE_FILE);

		Iterator<VariantContext> vcfIterator = SampleAssignmentVCFUtils.getVCFIterator(vcfReader, vcfSamples, !POLYMORPHIC_SNPS_ONLY, GQ_THRESHOLD,
				this.FRACTION_SAMPLES_PASSING, this.IGNORED_CHROMOSOMES, log);
		if (this.HET_SNPS_ONLY)
			vcfIterator = new HetSNPFilter(vcfIterator);

		SNPInfoCollection result = new SNPInfoCollection(vcfIterator, vcfReader.getFileHeader().getSequenceDictionary(), false, null, log);
		return (result);
	}

	private MultiCellDigitalAlleleCountsIterator getDACIterator(SNPUMIBasePileupIterator sbpi, boolean multiGenesPerRead, SNPInfoCollection snpInfo) {
		DigitalAlleleCountsGeneIteratorI dacIter = null;
		if (multiGenesPerRead) {
			dacIter = new DigitalAlleleCountsIterator(sbpi, BASE_QUALITY, snpInfo.getRefAllele(), snpInfo.getAltAllele());
		} else {
			dacIter = new DigitalAlleleCountsBestGeneIterator(sbpi, BASE_QUALITY, snpInfo.getRefAllele(), snpInfo.getAltAllele());
		}
		MultiCellDigitalAlleleCountsIterator multiIter = new MultiCellDigitalAlleleCountsIterator(dacIter);
		return (multiIter);
	}

	private PrintStream getPrintStreamOrNull(File outFile) {
		if (outFile == null)
			return null;
		PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(outFile));
		return (out);
	}

	/** Stock main method. */
	public static void main(final String[] args) {
		new GatherDigitalAlleleCounts().instanceMainWithExit(args);
	}
}
