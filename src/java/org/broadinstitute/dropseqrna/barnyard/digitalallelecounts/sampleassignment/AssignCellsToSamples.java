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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment;

import java.io.File;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import htsjdk.samtools.util.*;
import org.apache.commons.collections4.CollectionUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.BarcodeListRetrieval;
import org.broadinstitute.dropseqrna.barnyard.GeneFunctionCommandLineBase;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPInfoCollection;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPUMIBasePileup;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPUMIBasePileupIterator;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SortOrder;
import org.broadinstitute.dropseqrna.censusseq.JointIteratorCounter;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.AssertSequenceDictionaryIntersection;
import org.broadinstitute.dropseqrna.utils.Bases;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.IntervalTagComparator;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.SamHeaderUtil;
import org.broadinstitute.dropseqrna.utils.SamWriterSink;
import org.broadinstitute.dropseqrna.utils.VCFUtils;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.IgnoreGeneAnnotationTagger;
import org.broadinstitute.dropseqrna.utils.readiterators.PCRDuplicateFilteringIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.vcftools.SampleAssignmentVCFUtils;
import org.broadinstitute.dropseqrna.vcftools.filters.GenotypeGQFilter;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder.OutputType;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import picard.annotation.LocusFunction;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.BasicInputParser;
import picard.util.TabbedTextFileWithHeaderParser;
import picard.util.TabbedTextFileWithHeaderParser.Row;

@CommandLineProgramProperties(summary = "For a set of cells, and a VCF File," + "calculate the most likely sample assignment for each cell.  This is done by"
		+ "looking at each SNP/sample genotype, and calculating the likelihood of that result"
		+ "given the pileup of reads at that genomic location for each cell.  The log likelihoods"
		+ "across all genomic locations for each cell are calculated and summed.  "
		+ "Each cell then has a likelihood across all SNPs for each of the samples in the VCF file.  The most likely sample"
		+ "is assigned to the cell, and the likelihood ratio of the sample assignment to the next best sample assignment is "
		+ "calculated to give us a sense of how confident we are in the assignment.", oneLineSummary = "Assign cells in a BAM to the most likely sample in a VCF", programGroup = DropSeq.class)

public class AssignCellsToSamples extends GeneFunctionCommandLineBase {

	private static final Log log = Log.getInstance(AssignCellsToSamples.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze. This argument can accept wildcards, or a file with the suffix .bam_list that contains the locations of multiple BAM files", minElements = 1)
	public List<File> INPUT_BAM;

	@Argument(doc = "The input VCF file to analyze.")
	public File VCF;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output file of sample assignments. This supports zipped formats like gz and bz2.")
	public File OUTPUT;

	@Argument(doc = "Verbose output of every cell/SNP result. This supports zipped formats like gz and bz2.", optional = true)
	public File VERBOSE_OUTPUT;

	@Argument(doc = "Verbose output of every cell/SNP result for the BEST donor. This supports zipped formats like gz and bz2.", optional = true)
	public File VERBOSE_BEST_DONOR_OUTPUT = null;

	@Argument(doc = "Output a version of the VCF containing only the variants that have coverage in the input BAM.  All samples are retained.  This file can be used as the VCF input file on subsequent runs of the same data set to greatly speed up the run, if the same BAM is being analyzed with different conditions.  If you plan on calling doublets with DetectDoublets, you'll need this file.", optional = true)
	public File VCF_OUTPUT;

	@Argument(doc = "Output a version of the BAM containing only the reads that have coverage in the input BAM and VCF.", optional = true)
	public File BAM_OUTPUT;

	@Argument(doc = "The cell barcode tag.  If there are no reads with this tag, the program will assume that all reads belong to the same cell and process in single sample mode.")
	public String CELL_BARCODE_TAG = "XC";

	@Argument(doc = "The molecular barcode tag.")
	public String MOLECULAR_BARCODE_TAG = "XM";

	@Argument(doc = "The functional annotation for the read.  If set, extracts the functional annotation(s) [CODING/UTR/etc] at each SNP position and outputs in the verbose output.")
	public String FUNCTION_TAG = "XF";

	@Argument(doc = "The edit distance that molecular barcodes should be combined at within a gene/SNP.")
	public Integer EDIT_DISTANCE = 1;

	@Argument(doc = "The map quality of the read to be included.")
	public Integer READ_MQ = 10;

	@Argument(doc = "The minimum genotype quality for a variant.  Set this value to 0 to not filter by GQ scores if they are present, or to -1 to completely "
			+ "ignore GQ values if they are not set in the genotype info field.  If the GQ field is not set in the VCF header, this will be set to -1 by default.")
	public Integer GQ_THRESHOLD = 30;

	@Argument(doc = "Number of cells that you want to extract from the BAM. The program will picks the top <NUM_BARCODES> barcodes by read count.", mutex = {
			"CELL_BC_FILE" })
	public Integer NUM_BARCODES = null;

	@Argument(doc = "Override NUM_CORE_BARCODES and process reads that have the cell barcodes in this file instead.  The file has 1 column with no header.", mutex = {
			"NUM_BARCODES" })
	public File CELL_BC_FILE = null;

	@Argument(doc = "A file with a list of samples in the VCF to match up to cells in the BAM.  This subsets the VCF into a smaller data set containing only the samples listed. The file has 1 column with no header.", optional = true)
	public File SAMPLE_FILE = null;

	@Argument(doc = "This file provides an answer key for the input BAM and VCF.  This does not change any of the likelihood analysis results, but appends an additional column to the output that contains the known sample assignment for a cell."
			+ "This is useful when assessing how well the method works, but is completely optional.  The format of the file is 2 tab seperated columns with a header [cell 	sample].  The first column contains the cell barcodes from the BAM, "
			+ "the second the sample assignments from the VCF.", optional = true)
	public File ANSWER_KEY_FILE = null; 

	@Argument(doc = "A file that contains an estimate of how much ambient RNA is in each cell.  This is a fractional estimate between 0 and 1.  File is tab seperated, with 2 columns:"
			+ "cell_barcode and frac_contamination.  When supplied along side the ALLELE_FREQUENCY_ESTIMATE_FILE, this modifies the likelihood error rates to take into account how often"
			+ "the allele observed can be drawn from ambient RNA. We use cellbender remove background [https://github.com/broadinstitute/CellBender] to estimate the "
			+ "number of transcripts before and after ambient cleanup to define the fraction of transcripts that come from ambient RNA.", optional = true)
	public File CELL_CONTAMINATION_ESTIMATE_FILE = null;

	@Argument(doc = "A file that contains an estimate of the allele frequency expected for each SNP across donors. The best estimate of this will come from the fraction of reference and alternate allele"
			+ "UMIs that are observed at each snp site.  This report can be generated via GatherDigitalAlleleCounts.  This is a fractional estimate between 0 and 1.  File is tab seperated, with at least 3 columns:"
			+ "chromosome, position, maf_umi.  When supplied and CELL_CONTAMINATION_ESTIMATE_FILE is provided, this modifies the likelihood error rates to take into account how often "
			+ "the allele observed can be drawn from ambient RNA.", optional = true)
	public File ALLELE_FREQUENCY_ESTIMATE_FILE = null;

	private final String SNP_TAG = "YS";

	@Argument(doc = "Should monomorphic SNPs across the population be retained?")
	public boolean RETAIN_MONOMORPIC_SNPS = false;

	@Argument(doc = "At least <FRACTION_SAMPLES_PASSING> samples must have genotype scores >= GQ_THRESHOLD for the variant in the VCF to be included in the analysis.")
	public double FRACTION_SAMPLES_PASSING = 0.5;

	@Argument(doc = "A list of chromosomes to omit from the analysis.  The default is to omit the sex chromosomes.")
	public List<String> IGNORED_CHROMOSOMES = new ArrayList<>(Arrays.asList("X", "Y", "MT"));

	// to mimic the initial R implementation.
	@Argument(doc = "Instead of useing base qualities to determine error rate, use a fixed error rate instead. This is rounded to the nearest phread score internally.", optional = true)
	public Double FIXED_ERROR_RATE = null;

	@Argument(doc = "Caps the base error rate at a maximum probability so no SNP can be weighed more than this value.  For example, if this value was 0.01, "
			+ "then a base quality 30 value (normally an erro rate of 0.001) would become 0.01.  With the same threshold, a base with an error rate of 0.1 would be unaffected.", optional = true)
	public Double MAX_ERROR_RATE = null;

	// on a per-snp basis, generate a mixture of all the data across known samples to fill in likelihoods for the missing
	// samples.
	@Argument(doc = "on a per-snp basis, generate a mxiture of all the data across known samples to fill in likelihoods for the missing samples.")
	public Boolean ADD_MISSING_VALUES = true;

	private Boolean EMIT_POPULATION_LIKELIHOOD = true;

	@Argument(doc = "EXPERIMENTAL!!! Run the program in DNA Mode.  In this mode, reads should have a cell barcode, but will be missing gene annotations and UMIs.  All reads will be "
			+ "accepted as passing, and each read (or read pair) will be treated as a single UMI  If the data is PCR Duplicate marked, duplicate reads will be filtered. ")
	public Boolean DNA_MODE = false;

	@Argument(doc="Logging interval for SNP sequence pileup processing")
	public Integer SNP_LOG_RATE=1000;
	
	@Argument(doc="(Optional)  If set to a value, this program will quit early with an exception if this number of UMIs are observed without encountering a transcribed SNP - "
			+ "that is, a variant from the VCF file that is will be informative during donor assignment.  A value of -1 disables this test.  As some fraction of UMIs may not"
			+ "contain transcribed SNPs, the value for this should probably be set to a fairly large number, like 10% of the total UMIs in your experiment.  Mostly useful for debugging / or "
			+ "testing new data sets.")
	public Integer TRANSCRIBED_SNP_FAIL_FAST_THRESHOLD=-1;
	
	// This isn't fully implemented yet so is disabled.
	// @Argument(doc = "on a per-genotype basis, use the GQ score to weigh how much we believe the data given genotype.
	// Weights the likelihood (D|G) by the GQ score, which is converted from phred quality to a probability [0-1].")
	// public Boolean USE_GENOTYPE_QUALITY=false;

	// private Map <String,String> testMap = new HashMap<>();

	@Override
	protected int doWork() {
		// get the verbose output writer if not null, or return null.
		PrintStream verboseWriter = getVerboseWriter(this.VERBOSE_OUTPUT);

		final VCFFileReader vcfReader = new VCFFileReader(this.VCF, true);

		SamHeaderAndIterator headerAndIter = SamFileMergeUtil.mergeInputs(this.INPUT_BAM, false, SamReaderFactory.makeDefault());
		AssertSequenceDictionaryIntersection.assertIntersectionObjectVcf(headerAndIter.header, "BAM INPUT(S)", this.VCF, log);

		// disable GQ filter if it's not in the header.
		if (!VCFUtils.GQInHeader(vcfReader)) {
			this.GQ_THRESHOLD = -1;
			log.info("Genotype Quality [GQ] not found in header.  Disabling GQ_THRESHOLD parameter");
		}

		// subsetting the VCF by samples isn't generally needed for the standard way to run AssignCellsToSamples, but is super
		// expensive to run. We need to get around this.
		List<String> vcfSamples = SampleAssignmentVCFUtils.getVCFSamples(vcfReader, this.SAMPLE_FILE);

		// set up the VCF writer if needed.
		VariantContextWriter vcfWriter = getVcfWriter(vcfReader, vcfSamples);

		// populate the answer key if appropriate
		Map<String, String> answerKey = readAnswerKey();

		// set up the verbose output if needed
		writeVerboseOutputHeader(vcfSamples, answerKey, this.CELL_CONTAMINATION_ESTIMATE_FILE, this.ALLELE_FREQUENCY_ESTIMATE_FILE, verboseWriter);

		// get from the VCF!
		final SNPInfoCollection snpIntervals = SampleAssignmentVCFUtils.getSNPInfoCollection(this.VCF, vcfSamples, this.RETAIN_MONOMORPIC_SNPS, this.GQ_THRESHOLD,
				this.FRACTION_SAMPLES_PASSING, IGNORED_CHROMOSOMES, null, false, log);

		// a requested early exit if there are no SNPs.
		if (snpIntervals.isEmpty()) {
			log.error("No SNP intervals detected!  Check to see if your VCF filter thresholds are too restrictive!");
			return 1;
		}
			
		// this is the same method of filtering used to generate the snpIntervals. Don't ask for a progress logger, we'll use
		// the JointIteratorCounter below.
		// final PeekableIterator<VariantContext> vcfIterator = getVCFIterator (vcfReader, vcfSamples,
		// this.RETAIN_MONOMORPIC_SNPS, this.log);
		PeekableIterator<VariantContext> vcfIterator = SampleAssignmentVCFUtils.getVCFIterator(vcfReader, vcfSamples, this.RETAIN_MONOMORPIC_SNPS,
				this.GQ_THRESHOLD, this.FRACTION_SAMPLES_PASSING, IGNORED_CHROMOSOMES, log);

		// and the BAM iterator.
		List<String> barcodes = getCellBarcodes();
		PeekableIterator<List<SampleGenotypeProbabilities>> sampleGenotypeIterator = prepareIterator(snpIntervals, barcodes);

		// load up the contamination and MAF if available.
		Map<String, Double> contaminationMap = CellAssignmentUtils.getCellContamination(this.CELL_CONTAMINATION_ESTIMATE_FILE, barcodes);
		Map<Interval, Double> variantMinorAlleleFrequency = CellAssignmentUtils.getMinorAlleleFrequencyMap(this.ALLELE_FREQUENCY_ESTIMATE_FILE);

		// prepare the collection of likelihoods and other metrics
		CellCollectionSampleLikelihoodCollection likelihoodCollection = initializeLikelihoodCollection(vcfSamples, contaminationMap,
				variantMinorAlleleFrequency);

		// some counters for information on what's being iterated on. Every sample genotype result *should* match up to a VCF
		// entry
		// so we expect those numbers to be the same at the end.
		JointIteratorCounter counter = new JointIteratorCounter();

		// We need the sequence dictionary in the interval list for all sorting of chromosomes.
		SAMSequenceDictionary dict = snpIntervals.getIntervalList().getHeader().getSequenceDictionary();

		// TODO: should be able to plug VCFPileupJointIterator <T extends SNPIntervalRecordI>
		// Would have to wrap SampleGenotypeProbabilities in a collection that extends SNPIntervalRecordI for this to work,
		// yuck.

		// walk through both the VCF and the BAM at the same time.
		// the sampleGenotypeIterator gives all results for a SNP, so that's the results for every cell that had reads on the
		// SNP.
		while (vcfIterator.hasNext() && sampleGenotypeIterator.hasNext()) {
			List<SampleGenotypeProbabilities> pileUpsSingleSNP = sampleGenotypeIterator.peek();
			// only have to compare the first record.
			int cmp = compareRecords(pileUpsSingleSNP.get(0), vcfIterator.peek(), dict);
			if (cmp < 0) {
				sampleGenotypeIterator.next();
				counter.SAMPLE_GENE_ITER++;
			} else if (cmp > 0) {
				vcfIterator.next();
				counter.VCF_ITER++;
			} else if (cmp == 0) {
				counter.BOTH++;
				counter.SAMPLE_GENE_ITER++;
				counter.VCF_ITER++;
				if (counter.BOTH % this.SNP_LOG_RATE == 0)
					log.info("Processed [" + counter.BOTH + "] SNPs in BAM + VCF");

				// grab the next record and process it.
				VariantContext vc = vcfIterator.next();
				// write the VCF record out if needed.
				if (vcfWriter != null)
					vcfWriter.add(vc);
				pileUpsSingleSNP = sampleGenotypeIterator.next();
				likelihoodCollection = processSNP(vc, pileUpsSingleSNP, likelihoodCollection, ADD_MISSING_VALUES);
				if (verboseWriter != null) {
					writeVerboseOutput(vc, pileUpsSingleSNP, answerKey, vcfSamples, contaminationMap, variantMinorAlleleFrequency, verboseWriter);
				}

				if (counter.SAMPLE_GENE_ITER != counter.BOTH) {
					log.error("Expected counts don't match: SGICounter [" + counter.SAMPLE_GENE_ITER + "] both counter [" + counter.BOTH + "]");
					log.error(
							"VCF and BAM iteration have gone out of sync.  Iteration relies on the sequence dictonary order of the VCF file - make sure the body of the VCF is sorted in the order dictated by the VCF header.");
					vcfIterator.close();
					return 1;
				}
			}
		}
		log.info("Processed [" + counter.BOTH + "] SNPs in BAM + VCF");

		if (counter.SAMPLE_GENE_ITER != counter.BOTH)
			log.info("Expected counts don't match: SGICounter [" + counter.SAMPLE_GENE_ITER + "] both counter [" + counter.BOTH + "]");

		vcfIterator.close();
		vcfReader.close();
		if (vcfWriter != null)
			vcfWriter.close();
		sampleGenotypeIterator.close();

		if (this.VERBOSE_OUTPUT != null)
			verboseWriter.close();

		// write out results.
		ErrorCheckingPrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(this.OUTPUT));
		addMetaDataToHeader(out);
		writeBestLikelihoods(vcfSamples, likelihoodCollection, answerKey, out, this.EMIT_POPULATION_LIKELIHOOD);
		out.close();

		// filter out the VERBOSE results if required.
		if (this.VERBOSE_BEST_DONOR_OUTPUT != null) {
			log.info("Writing Verbose Best Donor Report");
			filterVerboseToBest(this.VERBOSE_OUTPUT, this.VERBOSE_BEST_DONOR_OUTPUT, likelihoodCollection);
		}
		log.info("Finished!");
		return 0;
	}

	// helper method to deal with parsing input file if it exists.
	private CellCollectionSampleLikelihoodCollection initializeLikelihoodCollection(List<String> vcfSamples, final Map<String, Double> contaminationMap,
			Map<Interval, Double> variantMinorAlleleFrequency) {
		if (contaminationMap == null)
			return (new CellCollectionSampleLikelihoodCollection(FIXED_ERROR_RATE, this.MAX_ERROR_RATE, vcfSamples));
		return new CellCollectionSampleLikelihoodCollection(contaminationMap, variantMinorAlleleFrequency, vcfSamples);
	}

	/**
	 * Takes the
	 * 
	 * @param verboseOutput
	 * @param bestDonorVerboseOutput
	 * @param likelihoodCollection
	 */
	static void filterVerboseToBest(final File verboseOutput, final File bestDonorVerboseOutput,
			CellCollectionSampleLikelihoodCollection likelihoodCollection) {
		Map<String, String> bestDonorPerCell = new HashMap<String, String>();
		// get the best donor per cell.
		likelihoodCollection.getCellBarcodes().stream()
				.forEach(x -> bestDonorPerCell.put(x, likelihoodCollection.getLikelihoodCollection(x).getBestSampleAssignment().getSample()));

		// get the list of all samples
		// Collection<String> donorList = likelihoodCollection.getSamples();

		TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(verboseOutput);
		// Set<String> colNames = parser.getColumnNames();

		CloseableIterator<Row> rowIter = parser.iterator();
		log.info("Processing verbose donor assignment file to generate best donor output");
		ProgressLogger pl = new ProgressLogger(log, 1000000, "Processed");

		// prepare the writer.
		PrintStream writer = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(bestDonorVerboseOutput));

		// fixed fields we want to extract.
		List<String> fixedFields = new ArrayList<>();
		fixedFields.addAll(Arrays.asList("chr", "pos", "gene", "function", "cell"));
		for (Bases b : Bases.values())
			fixedFields.add(b.toString());
		fixedFields.addAll(Arrays.asList("num_umi", "REF", "ALT"));

		// the output header
		List<String> outputHeader = new ArrayList<>(fixedFields);

		boolean outputPopAverage = false;
		boolean outputMaf = false;
		if (parser.getColumnNames().contains("population_average_likelihood"))
			outputPopAverage = true;
		if (parser.getColumnNames().contains("minor_allele_freq"))
			outputMaf = true;

		outputHeader.addAll(Arrays.asList("best_donor", "genotype", "log10_likelihood"));
		if (outputPopAverage)
			outputHeader.add("population_average_likelihood");
		if (outputMaf)
			outputHeader.add("minor_allele_freq");

		writer.println(StringUtils.join(outputHeader, "\t"));

		// write the body.
		while (rowIter.hasNext()) {
			Row r = rowIter.next();
			pl.record(r.getField("chr"), Integer.parseInt(r.getField("pos")));

			List<String> line = new ArrayList<>();
			// get all the fixed fields into the output.
			for (String field : fixedFields) {
				line.add(r.getField(field));
			}
			// get the best donor for this cell and output the likelihood for their genotype
			String cell = r.getField("cell");
			String bestDonor = bestDonorPerCell.get(cell);
			String genotype = r.getField(bestDonor);
			// This is slightly lame, but I'm currently encoding genotypes as "null" in the verbose output, and
			// that's getting parsed as an actual null, so need to catch it.
			if (genotype == null)
				genotype = "NA";
			line.add(bestDonor);
			line.add(genotype);

			switch (genotype) {
			case "ref":
				line.add(r.getField("ref_loglikelihood"));
				break;
			case "het":
				line.add(r.getField("het_loglikelihood"));
				break;
			case "alt":
				line.add(r.getField("alt_loglikelihood"));
				break;
			default:
				line.add("NA");
			}
			if (outputPopAverage)
				line.add(r.getField("population_average_likelihood"));
			if (outputMaf)
				line.add(r.getField("minor_allele_freq"));
			writer.println(StringUtils.join(line, "\t"));
		}
		CloserUtil.close(parser);
		writer.close();
	}

	static void writeBestLikelihoods(final List<String> samples, final CellCollectionSampleLikelihoodCollection likelihoodCollection,
			final Map<String, String> answerKey, final ErrorCheckingPrintStream out, boolean emitPopulationLikelihood) {
		// write header
		writeBestLikelihoodsHeader(samples, answerKey, out, emitPopulationLikelihood);

		// write body
		DecimalFormat likeFormat = new DecimalFormat("0.####E0");

		Map<String, Double> fdr = likelihoodCollection.getFDRCorrectedPvalues();

		// int counter=0;
		for (String cellBarcode : likelihoodCollection.getCellBarcodes()) {

			CellSampleLikelihoodCollection cslc = likelihoodCollection.getLikelihoodCollection(cellBarcode);
			List<String> line = new ArrayList<>();
			line.add(cellBarcode);
			// the number of snp and umi observations
			line.add(Integer.toString(cslc.getNumSNPs()));
			line.add(Integer.toString(cslc.getNumUMIs()));

			// write best sample info.
			BestSampleAssignmentForCell b = cslc.getBestSampleAssignment();
			line.add(likeFormat.format(b.getLogLikelihoodRatio()));
			line.add(likeFormat.format(b.getOneMinusPvalue()));
			line.add(likeFormat.format(fdr.get(cellBarcode)));
			line.add(likeFormat.format(b.getBestLoglikelihood()));
			line.add(b.getSample());
			line.add(likeFormat.format(cslc.getMedianLikelihood()));
			if (emitPopulationLikelihood)
				line.add(likeFormat.format(cslc.getGlobalPenaltyScore()));
			if (answerKey != null) {
				String answer = answerKey.get(cellBarcode);
				if (answer == null)
					line.add("NA");
				else
					line.add(answer);
			}

			// add each sample's likelihood.

			for (String sample : samples)
				if (cslc.isSampleInCollection(sample)) {
					Double l = cslc.getLoglikelihood(sample);
					line.add(likeFormat.format(l));
				} else
					line.add("NA");

			String h = StringUtils.join(line, "\t");
			out.println(h);
		}
	}

	private void addMetaDataToHeader(final PrintStream out) {

		List<String> paths = this.INPUT_BAM.stream().map(x -> x.getAbsolutePath()).collect(Collectors.toList());

		String bamList = StringUtils.join(paths, ",");
		List<String> header = new ArrayList<>(Arrays.asList("#INPUT_BAM=" + bamList, "INPUT_VCF=" + this.VCF.getAbsolutePath(),
				"DONOR_FILE=" + this.SAMPLE_FILE, "CELL_BC_FILE=" + CELL_BC_FILE, "GQ_THRESHOLD=" + Integer.toString(this.GQ_THRESHOLD),
				"FRACTION_SAMPLES_PASSING=" + Double.toString(this.FRACTION_SAMPLES_PASSING), "READ_MQ=" + Integer.toString(this.READ_MQ),
				"FIXED_ERROR_RATE=" + CellAssignmentUtils.convertNullToString(this.FIXED_ERROR_RATE),
				"MAX_ERROR_RATE=" + CellAssignmentUtils.convertNullToString(this.MAX_ERROR_RATE), "LOCUS_FUNCTION=" + this.LOCUS_FUNCTION_LIST.toString()));

		if (this.CELL_CONTAMINATION_ESTIMATE_FILE != null && this.ALLELE_FREQUENCY_ESTIMATE_FILE != null) {
			header.add("CELL_CONTAMINATION_ESTIMATE_FILE=" + CELL_CONTAMINATION_ESTIMATE_FILE.getAbsolutePath());
			header.add("ALLELE_FREQUENCY_ESTIMATE_FILE=" + ALLELE_FREQUENCY_ESTIMATE_FILE.getAbsolutePath());
		}

		String h = StringUtils.join(header, "\t");
		out.println(h);
	}

	static void writeBestLikelihoodsHeader(final List<String> samples, final Map<String, String> answerKey, final PrintStream out,
			boolean emitPopulationLikelihood) {
		List<String> line = new ArrayList<>();
		line.add("cell");
		line.add("num_snps");
		line.add("num_umis");
		line.add("ratio");
		line.add("pvalue");
		line.add("FDR_pvalue");
		line.add("bestLikelihood");
		line.add("bestSample");
		line.add("median_likelihood");
		if (emitPopulationLikelihood)
			line.add("population_average_likelihood");
		if (answerKey != null)
			line.add("trueSample");

		line.addAll(samples);
		String h = StringUtils.join(line, "\t");
		out.println(h);
	}

	CellCollectionSampleLikelihoodCollection processSNP(final VariantContext vc, final List<SampleGenotypeProbabilities> sgpList,
			final CellCollectionSampleLikelihoodCollection likelihoodCollection, final boolean addMissingValues) {
		Set<String> sampleNames = vc.getSampleNames();
		return processSNP(vc, sgpList, sampleNames, likelihoodCollection, addMissingValues);
	}

	/**
	 * @param b byte value of A,C, G or T
	 * @return very compact representation of A, C, G or T (2 bits)
	 */
	private static int getCompactBaseRepresentation(byte b) {
		switch (b) {
			case SequenceUtil.A:
				return 0;
			case SequenceUtil.C:
				return 1;
			case SequenceUtil.G:
				return 2;
			case SequenceUtil.T:
				return 3;
			default:
				throw new IllegalArgumentException(String.format("Unexpected base %d", b));
		}
	}

	/**
	 *
	 * @param a1 byte value of A,C, G or T
	 * @param a2 byte value of A,C, G or T
	 * @return very compact representation of this ordered pair of alleles (4 bits)
	 */
	private static int getAllelesToGroupBy(final byte a1, final byte a2) {
		return (getCompactBaseRepresentation(a1) << 2) | getCompactBaseRepresentation(a2);
	}

	// Return pair of alleles in compact form for Collectors.groupingBy
	private static int getAllelesToGroupBy(final Genotype g) {
		return getAllelesToGroupBy(g.getAllele(0).getBases()[0], g.getAllele(1).getBases()[0]);
	}

	// We know that T,T is the largest value in the compact rep.
	private static final int NumCompactAlleleRepresentations =
			getAllelesToGroupBy(SequenceUtil.T, SequenceUtil.T) + 1;

	// Recycle to avoid GC
	// The index into these arrays is produced by getAllelesToGroupBy()
	private final Genotype[] prototypeGenotypeForAllelePair = new Genotype[NumCompactAlleleRepresentations];
	// Once a list is created, it is not destroyed but rather cleared for the next SNP.
	// However, if prototypeGenotypeForAllelePair[i] is null, then samplesForAllelePair[i] may contain sample names
	// from a previous iteration that are no longer relevant.
	private final List<String>[] samplesForAllelePair = new List[NumCompactAlleleRepresentations];

	/**
	 * This updates the likelihood all cells that have one or more UMIs at this variant.
	 * 
	 * @param vc
	 * @param sgpList
	 * @param sampleNames
	 * @param likelihoodCollection
	 * @param addMissingValues
	 * @return
	 */
	CellCollectionSampleLikelihoodCollection processSNP(final VariantContext vc, final List<SampleGenotypeProbabilities> sgpList, final Set<String> sampleNames,
			CellCollectionSampleLikelihoodCollection likelihoodCollection, final boolean addMissingValues) {
		GenotypesContext gc = vc.getGenotypes(sampleNames);
		Iterator<Genotype> gi = gc.iterator();
		GenotypeGQFilter gf = new GenotypeGQFilter(gi, this.GQ_THRESHOLD);
		char refAlleleBase = StringUtil.byteToChar(vc.getReference().getBases()[0]);

		// Primitive implementation of grouping together the samples that have the same pair of alleles for this SNP.
		// This way has many fewer objects created than groupBy.

		// This array holds the first Genotype object encountered for a pair of alleles.  The index into the
		// array is produced by getAllelesToGroupBy()
		Arrays.fill(prototypeGenotypeForAllelePair, null);
		// Don't bother clearing samplesForAllelePair globally.  Just clear the entries for allele pairs
		// encountered for this SNP.

		// Capture one Genotype object for each distinct pair of alleles, and get the list of
		// sample names for each distinct pair of alleles.
		for (final Genotype g: gf) {
			final int compactAlleleRep = getAllelesToGroupBy(g);
			if (prototypeGenotypeForAllelePair[compactAlleleRep] == null) {
				prototypeGenotypeForAllelePair[compactAlleleRep] = g;
				if (samplesForAllelePair[compactAlleleRep] == null) {
					samplesForAllelePair[compactAlleleRep] = new ArrayList<>(sampleNames.size());
				} else {
					samplesForAllelePair[compactAlleleRep].clear();
				}
			}
			samplesForAllelePair[compactAlleleRep].add(g.getSampleName());
		}

		// batch update all donors that have the same genotype.
		// map allele pair to list of Genotype.
		// then batch update likelihoods on the list of sample names.
		// this would be a bad plan if we used genotype quality, where each sample could have a different genotype quality
		// score.
		for (int i = 0; i < NumCompactAlleleRepresentations; ++i) {
			if (prototypeGenotypeForAllelePair[i] == null) {
				continue;
			}
			final Genotype g = prototypeGenotypeForAllelePair[i];
			byte[] a1Bases = g.getAllele(0).getBases();
			byte[] a2Bases = g.getAllele(1).getBases();
			char a1 = StringUtil.byteToChar(a1Bases[0]);
			char a2 = StringUtil.byteToChar(a2Bases[0]);
			final List<String> sampleNamesForAlleles = samplesForAllelePair[i];
			// this isn't fully implemented yet so is disabled.
			/*
			 * Double genotypeProbability=null; if (USE_GENOTYPE_QUALITY) { Integer genotypeQuality = g.getGQ();
			 * genotypeProbability=LikelihoodUtils.getInstance().phredScoreToErrorProbability(genotypeQuality.byteValue()); }
			 */
			// for this sample genotype, update all the cells where we have data.
			for (SampleGenotypeProbabilities p : sgpList) {
				likelihoodCollection.updatelikelihoods(p, sampleNamesForAlleles, a1, a2, null, refAlleleBase);
			}

		}

		CloserUtil.close(gf);

		// update number of SNPs and umis observed once per SNP in the VCF.
		for (SampleGenotypeProbabilities p : sgpList)
			likelihoodCollection.updateNumObservations(p);

		// calculate the global penalty score. Optionally apply this score to donor likelihoods for a cell/snp where the donor
		// does not have a genotype set.
		likelihoodCollection = calculateMissingLikelihoodOptionallyApply(vc, sgpList, sampleNames, likelihoodCollection, addMissingValues);
		return likelihoodCollection;

	}

	/**
	 * For a SNP, runs this process against all cells:
	 * 
	 * Calculates the population-wide likelihood score, computed as the weighted average likelihood of the hom/het/var
	 * likelihoods weighted by all the donors at this SNP that have valid genotypes This population wide penalty score is
	 * incremented for every cell Optionally, this score is used as the likelihood score for all donors that do not have a
	 * valid genotype at this SNP.
	 * 
	 * @param vc                   The current SNP from the VCF
	 * @param sgpList              The list of SNP pileups objects - one per cell.
	 * @param sampleNames          The original set of samples to process from the VCF.
	 * @param likelihoodCollection The set of likelihood scores that are iteratively computed. This object will be updated
	 *                             and returned.
	 * @param applyToMissingDonors Set to true to assign the population wide likelihood score to all donors that don't have
	 *                             a genotype at this SNP.
	 * @return The updated likelihood scores and other tracking metrics across snps / donors / cells.s
	 */
	CellCollectionSampleLikelihoodCollection calculateMissingLikelihoodOptionallyApply(final VariantContext vc, final List<SampleGenotypeProbabilities> sgpList,
			final Set<String> sampleNames, final CellCollectionSampleLikelihoodCollection likelihoodCollection, boolean applyToMissingDonors) {

		// this class simplifies out a lot of code that was here.
		GenotypeCountMetrics gm = new GenotypeCountMetrics(vc, sampleNames, this.GQ_THRESHOLD);
		Double genotypeProbability = null;

		/*
		 * if (this.USE_GENOTYPE_QUALITY) { // the mean of the genotype qualities? if (meanGenotypeQuality.getN()==0)
		 * log.info("STOP"); double meanGQ=gm.getMeanGenotypeQuality();
		 * genotypeProbability=LikelihoodUtils.getInstance().phredScoreToErrorProbability(new Double(meanGQ).byteValue()); }
		 */

		// for each cell, update the samples where we have no genotype assigned for this variant.
		for (SampleGenotypeProbabilities p : sgpList) {
			double globalPenaltyLogLikelihood = likelihoodCollection.getMissingLikelihood(p, gm.getRefAllele(), gm.getAltAllele(), gm.getRefCount(),
					gm.getHetCount(), gm.getAltCount(), genotypeProbability);
			likelihoodCollection.incrementGlobalPenaltyAllDonors(p, globalPenaltyLogLikelihood);
			if (applyToMissingDonors)
				likelihoodCollection.applyMissingLikelihoodToDonors(p, globalPenaltyLogLikelihood, gm.getSamplesNotObserved());
		}
		return (likelihoodCollection);
	}

	private void writeVerboseOutputHeader(final List<String> sampleNames, final Map<String, String> answerKey, File contaminationMap,
			File variantMinorAlleleFrequency, final PrintStream out) {
		if (out == null)
			return;
		addMetaDataToHeader(out);
		List<String> line = new ArrayList<>();
		line.add("chr");
		line.add("pos");
		line.add("gene");
		if (this.FUNCTION_TAG != null)
			line.add("function");
		line.add("cell");
		for (Bases b : Bases.values())
			line.add(b.toString());
		line.add("num_umi");
		line.add("REF");
		line.add("ALT");
		line.add("ref_loglikelihood");
		line.add("het_loglikelihood");
		line.add("alt_loglikelihood");
		line.add("population_average_likelihood");
		line.add("minor_allele_freq");
		if (contaminationMap != null)
			line.add("contamination");
		if (variantMinorAlleleFrequency != null)
			line.add("umi_maf");

		if (answerKey != null)
			line.add("knownSample");

		for (String sample : sampleNames)
			line.add(sample);

		String h = StringUtils.join(line, "\t");
		out.println(h);
	}

	private String getLocusFunctionString(final SampleGenotypeProbabilities p) {
		Set<LocusFunction> lfSet = p.getLocusFunctions();
		if (lfSet == null || lfSet.size() == 0)
			return "";

		StringBuilder b = new StringBuilder();

		Iterator<LocusFunction> iter = lfSet.iterator();
		b.append(iter.next());

		while (iter.hasNext()) {
			b.append(",");
			b.append(iter.next());
		}
		return b.toString();

	}

	private void writeVerboseOutput(final VariantContext vc, final List<SampleGenotypeProbabilities> sgpList, final Map<String, String> answerKey,
			final List<String> sampleNames, Map<String, Double> contaminationMap, Map<Interval, Double> variantMinorAlleleFrequency, final PrintStream out) {
		if (out == null)
			return;
		DecimalFormat likeFormat = new DecimalFormat("0.######################");
		DecimalFormat mafFormat = new DecimalFormat("0.###");

		// This is expensive, but I'm not sure of a better way to get at the global penalty score per variant.
		// TODO: rework the code so this doesn't need to be computed once during aggregation and once here.
		CellCollectionSampleLikelihoodCollection likelihoodCollection = new CellCollectionSampleLikelihoodCollection(FIXED_ERROR_RATE, this.MAX_ERROR_RATE,
				sampleNames);
		likelihoodCollection = processSNP(vc, sgpList, likelihoodCollection, ADD_MISSING_VALUES);

		// get genotype information
		GenotypeCountMetrics gm = new GenotypeCountMetrics(vc, sampleNames, this.GQ_THRESHOLD);
		double maf = gm.getMinorAlleleFrequency();

		// each cell is one line for the same SNP.
		for (SampleGenotypeProbabilities p : sgpList) {
			List<String> line = new ArrayList<>();

			Double contamination = CellAssignmentUtils.getNullableValue(contaminationMap, p.getCell());
			Double umiMaf = CellAssignmentUtils.getNullableValue(variantMinorAlleleFrequency, p.getSNPInterval());

			line.add(p.getSNPInterval().getContig());
			line.add(Integer.toString(p.getSNPInterval().getStart()));
			
			/*
			for (SNPUMIBasePileup pileup: p.getBackingPileups()) {
				log.info(pileup);
			}
			*/
			
			String gene = p.getBackingPileups().iterator().next().getGene();
			line.add(gene);

			if (this.FUNCTION_TAG != null)
				line.add(getLocusFunctionString(p));

			line.add(p.getCell());

			// iterate over bases and get #umis
			ObjectCounter<Character> umiCount = p.getUMIBaseCounts();

			for (Bases b : Bases.values()) {
				int n = umiCount.getCountForKey(b.getBase());
				line.add(Integer.toString(n));
			}
			line.add(Integer.toString(umiCount.getTotalCount()));
			// these work only because we only allow single base reference/alternate alleles
			char refBase = StringUtil.byteToChar(vc.getReference().getBases()[0]);
			char altBase = CellAssignmentUtils.getAltAllele(vc);

			line.add(String.valueOf(refBase));
			line.add(String.valueOf(altBase));

			// TODO is there a way to plug in if the genotype quality is used? It's hardcoded off right now, but the feature is
			// disabled.
			double likelihoodRef = p.getLogLikelihood(refBase, refBase, this.FIXED_ERROR_RATE, null, this.MAX_ERROR_RATE, refBase, umiMaf, contamination);
			double likelihoodHet = p.getLogLikelihood(refBase, altBase, this.FIXED_ERROR_RATE, null, this.MAX_ERROR_RATE, refBase, umiMaf, contamination);
			double likelihoodAlt = p.getLogLikelihood(altBase, altBase, this.FIXED_ERROR_RATE, null, this.MAX_ERROR_RATE, refBase, umiMaf, contamination);

			line.add(likeFormat.format(likelihoodRef));
			line.add(likeFormat.format(likelihoodHet));
			line.add(likeFormat.format(likelihoodAlt));

			double populationPenalty = likelihoodCollection.getLikelihoodCollection(p.getCell()).getGlobalPenaltyScore();
			line.add(likeFormat.format(populationPenalty));
			line.add(mafFormat.format(maf));

			if (contaminationMap != null)
				line.add(mafFormat.format(contamination));
			if (variantMinorAlleleFrequency != null)
				line.add(mafFormat.format(umiMaf));

			if (answerKey != null) {
				String sample = answerKey.get(p.getCell());
				line.add(sample);
			}

			List<String> genoStrings = sampleNames.stream().map(x -> gm.getGenotypeString(x)).collect(Collectors.toList());
			line.addAll(genoStrings);

			String h = StringUtils.join(line, "\t");
			out.println(h);
		}

	}

	/**
	 * Compare the genotype probabilities location to the variant context probability. return 0 < if the sgp is before than
	 * the VC return 0 if the GP = VC return > 0 if the VC is before the sgp This is a comparison of the interval of the
	 * SampleGenotypeProbabilities to the VariantContext: sgp.compareTo(vc)
	 * 
	 * @param sgp
	 * @param vc
	 * @return
	 */
	int compareRecords(final SampleGenotypeProbabilities sgp, final VariantContext vc, final SAMSequenceDictionary dict) {
		Interval sgpI = sgp.getSNPInterval();
		Interval sgpVC = new Interval(vc.getContig(), vc.getStart(), vc.getEnd());
		// log.info("SNP Interval" + sgpI.toString() + " VCF interval " + sgpVC.toString());
		int cmp = IntervalTagComparator.compare(sgpI, sgpVC, dict);
		return cmp;
	}

	// done before, pretty well tested elsewhere boilerplate.
	public PeekableIterator<List<SampleGenotypeProbabilities>> prepareIterator(final SNPInfoCollection snpIntervals, List<String> barcodes) {
		// generate the reader with the desired options so we can enable EAGERLY_DECODE.
		SamReaderFactory factory = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE);
		SamHeaderAndIterator headerAndIter = SamFileMergeUtil.mergeInputs(this.INPUT_BAM, false, factory);

		// override the normal gene annotations with new ones before any other operations.
		// filter out PCR duplicates.
		if (this.DNA_MODE) {
			// Replace UMI tags with read names, set strand tag to read tag, set gene function to an accepted on. Overwrites ALL
			// reads tags.
			IgnoreGeneAnnotationTagger tagger = new IgnoreGeneAnnotationTagger(headerAndIter.iterator, this.GENE_NAME_TAG, this.GENE_STRAND_TAG,
					this.GENE_FUNCTION_TAG, this.LOCUS_FUNCTION_LIST, false, this.MOLECULAR_BARCODE_TAG, true);
			headerAndIter = new SamHeaderAndIterator(headerAndIter.header, (CloseableIterator<SAMRecord>) tagger.iterator());
			final PCRDuplicateFilteringIterator pcrDuplicateFilteringIterator = new PCRDuplicateFilteringIterator(headerAndIter.iterator);
			headerAndIter = new SamHeaderAndIterator(headerAndIter.header, pcrDuplicateFilteringIterator);
		}
		
		// this explicitly filters to the best SNP seen on a read if there are multiple SNPs touched by a read
		Map<Interval, Double> genotypeQuality = snpIntervals.getAverageGQ();		
		SNPUMIBasePileupIterator sbpi = new SNPUMIBasePileupIterator(headerAndIter, snpIntervals.getIntervalList(), GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG,
				LOCUS_FUNCTION_LIST, STRAND_STRATEGY, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, this.SNP_TAG, this.FUNCTION_TAG, this.READ_MQ, false,
				barcodes, genotypeQuality, SortOrder.SNP_CELL, this.TRANSCRIBED_SNP_FAIL_FAST_THRESHOLD);
		
		
		if (this.BAM_OUTPUT != null) {
			SAMFileHeader header = headerAndIter.header;
			SamHeaderUtil.addPgRecord(header, this);
			header.addComment("Reads that overlap SNPs in VCF [" + this.VCF + "]");
			SAMFileWriter writer = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, false, this.BAM_OUTPUT);
			final SamWriterSink informativeRecSink = new SamWriterSink(writer);
			sbpi.addReadSink(informativeRecSink);
		}

		final SAMSequenceDictionary dict = snpIntervals.getIntervalList().getHeader().getSequenceDictionary();

		// gets a SampleGenotypeProbabilities for each cell.
		SampleGenotypeProbabilitiesIterator result = new SampleGenotypeProbabilitiesIterator(sbpi, dict, this.EDIT_DISTANCE, SortOrder.SNP_CELL);

		// clusters SampleGenotypeProbabilities objects across all cells for a SNP
		GroupingIterator<SampleGenotypeProbabilities> groupingIterator = new GroupingIterator<>(result, new Comparator<SampleGenotypeProbabilities>() {
			@Override
			public int compare(final SampleGenotypeProbabilities o1, final SampleGenotypeProbabilities o2) {
				int cmp = IntervalTagComparator.compare(o1.getSNPInterval(), o2.getSNPInterval(), dict);
				return cmp;
			}
		});

		PeekableIterator<List<SampleGenotypeProbabilities>> peekableIter = new PeekableIterator<>(groupingIterator);
		return (peekableIter);
	}

	// done before, pretty well tested elsewhere boilerplate.
	public List<String> getCellBarcodes() {

		if (this.CELL_BC_FILE != null) {
			IOUtil.assertFileIsReadable(this.CELL_BC_FILE);
			List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(this.CELL_BC_FILE);
			log.info("Found " + cellBarcodes.size() + " cell barcodes in file");
			return (cellBarcodes);
		}
		log.info("Gathering barcodes for the top [" + this.NUM_BARCODES + "] cells");
		return new BarcodeListRetrieval().getListCellBarcodesByReadCount(this.INPUT_BAM, this.CELL_BARCODE_TAG, this.READ_MQ, null, this.NUM_BARCODES);
	}

	public Map<String, String> readAnswerKey() {
		if (this.ANSWER_KEY_FILE == null)
			return null;
		IOUtil.assertFileIsReadable(this.ANSWER_KEY_FILE);

		Map<String, String> result = new HashMap<>();

		BasicInputParser parser = new BasicInputParser(false, 2, this.ANSWER_KEY_FILE);
		if (parser.hasNext()) {
			String[] line = parser.next();
			if (!line[0].equals("cell") || !line[1].equals("sample")) {
				parser.close();
				throw new IllegalArgumentException("Wrong header, is this the answer key file you're looking for? " + parser.getCurrentLine());
			}
		}
		while (parser.hasNext()) {
			String[] line = parser.next();
			if (line.length != 2) {
				parser.close();
				throw new IllegalArgumentException("Expected answer key file to have 2 columns per line, but this line did not " + parser.getCurrentLine());
			}
			result.put(line[0], line[1]);
		}
		parser.close();
		return (result);
	}

	/**
	 * Returns a writer, or null if the out file is null.
	 * 
	 * @param vcfReader
	 * @param out
	 * @return
	 */
	VariantContextWriter getVCFWriter(final VCFFileReader vcfReader, final File out) {
		if (out == null)
			return null;
		VariantContextWriterBuilder builder = new VariantContextWriterBuilder().setReferenceDictionary(vcfReader.getFileHeader().getSequenceDictionary())
				.setOption(Options.INDEX_ON_THE_FLY).setBuffer(8192).setOutputFile(out);

		if (out.getName().endsWith(".gz"))
			builder.setOutputFileType(OutputType.BLOCK_COMPRESSED_VCF);
		else
			builder.setOutputFileType(OutputType.VCF);

		VariantContextWriter writer = builder.build();

		return writer;
	}

	/**
	 * If there's a file with a list of samples to use, parse that and use that subset of samples Otherwise, get the list of
	 * samples in the VCF header.
	 * 
	 * @param vcfReader
	 * @return
	 */
	public List<String> getVCFSamples(final VCFFileReader vcfReader) {
		if (this.SAMPLE_FILE == null)
			return vcfReader.getFileHeader().getSampleNamesInOrder();
		IOUtil.assertFileIsReadable(this.SAMPLE_FILE);
		return ParseBarcodeFile.readCellBarcodeFile(SAMPLE_FILE);
	}

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new AssignCellsToSamples().instanceMain(args));
	}

	/**
	 * A helper class to gather statistics about one line of a VCF file
	 * 
	 * @author nemesh
	 *
	 */
	private class GenotypeCountMetrics {
		private final int genotypeQualityThreshold;
		private final Set<String> sampleNamesRequested;
		private double meanGenotypeQuality;
		private char refAllele;
		private char altAllele;
		private ObjectCounter<GenotypeType> genotypeCounts;
		private Map<String, GenotypeType> genotypes;

		public GenotypeCountMetrics(final VariantContext vc, final Collection<String> sampleNames, int genotypeQualityThreshold) {
			this.genotypeQualityThreshold = genotypeQualityThreshold;
			this.sampleNamesRequested = new HashSet<String>(sampleNames);
			this.genotypes = new HashMap<>();

			GenotypesContext gc = vc.getGenotypes(sampleNamesRequested);
			Iterator<Genotype> gi = gc.iterator();
			GenotypeGQFilter gf = new GenotypeGQFilter(gi, genotypeQualityThreshold);

			Mean meanGenotypeQuality = new Mean();

			this.genotypeCounts = new ObjectCounter<>();

			for (Genotype g : gf) {
				genotypeCounts.increment(g.getType());
				genotypes.put(g.getSampleName(), g.getType());
				meanGenotypeQuality.increment(g.getGQ());
			}

			gf.close();
			this.meanGenotypeQuality = meanGenotypeQuality.getResult();
			this.refAllele = StringUtil.byteToChar(vc.getReference().getBases()[0]);
			this.altAllele = CellAssignmentUtils.getAltAllele(vc);
		}

		public GenotypeType getGenotype(String sample) {
			return (this.genotypes.get(sample));
		}

		public Collection<String> getSamplesNotObserved() {
			return (CollectionUtils.subtract(this.sampleNamesRequested, this.getSampleNamesObserved()));
		}

		public int getRefCount() {
			return this.genotypeCounts.getCountForKey(GenotypeType.HOM_REF);
		}

		public int getHetCount() {
			return genotypeCounts.getCountForKey(GenotypeType.HET);
		}

		public int getAltCount() {
			return genotypeCounts.getCountForKey(GenotypeType.HOM_VAR);
		}

		public double getMinorAlleleFrequency() {
			// stupid java int vs double.
			double ref = getRefCount();
			double het = getHetCount();
			double alt = getAltCount();

			double minorAlleleFreq = ((alt * 2) + het) / (ref * 2 + het * 2 * alt * 2);
			return (minorAlleleFreq);
		}

		@SuppressWarnings("unused")
		public int getGenotypeQualityThreshold() {
			return genotypeQualityThreshold;
		}

		@SuppressWarnings("unused")
		public Set<String> getSampleNamesRequested() {
			return sampleNamesRequested;
		}

		public Collection<String> getSampleNamesObserved() {
			return this.genotypes.keySet();
		}

		@SuppressWarnings("unused")
		public double getMeanGenotypeQuality() {
			return meanGenotypeQuality;
		}

		public char getRefAllele() {
			return refAllele;
		}

		public char getAltAllele() {
			return altAllele;
		}

		public String getGenotypeString(String sample) {
			GenotypeType gt = this.getGenotype(sample);
			String result = "null";
			if (gt == null) {
				return (result);
			} else {
				switch (gt) {
				case HOM_REF:
					result = "ref";
					break;
				case HET:
					result = "het";
					break;
				case HOM_VAR:
					result = "alt";
					break;
				default:
					result = "NA";
				}
			}
			return (result);
		}
	}

	private VariantContextWriter getVcfWriter(final VCFFileReader vcfReader, List<String> vcfSamples) {
		// set up the VCF writer if needed.
		VariantContextWriter vcfWriter = null;
		if (this.VCF_OUTPUT != null) {
			IOUtil.assertFileIsWritable(this.VCF_OUTPUT);
			vcfWriter = SampleAssignmentVCFUtils.getVCFWriter(vcfReader, this.VCF_OUTPUT);
			VCFHeader header = vcfReader.getFileHeader();
			VCFHeader h2 = new VCFHeader(header.getMetaDataInInputOrder(), vcfSamples);
			vcfWriter.writeHeader(h2);
		}
		return (vcfWriter);
	}

	private PrintStream getVerboseWriter(File output) {
		PrintStream verboseWriter = null;
		if (output != null) {
			IOUtil.assertFileIsWritable(output);
			verboseWriter = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(output));
		}
		return (verboseWriter);
	}

	@Override
	protected String[] customCommandLineValidation() {

		final ArrayList<String> list = new ArrayList<>(1);

		IOUtil.assertFileIsReadable(this.VCF);
		this.INPUT_BAM = FileListParsingUtils.expandFileList(INPUT_BAM);
		IOUtil.assertFileIsWritable(this.OUTPUT);
		if (CELL_BC_FILE != null)
			IOUtil.assertFileIsReadable(this.CELL_BC_FILE);
		if (this.BAM_OUTPUT != null)
			IOUtil.assertFileIsWritable(this.BAM_OUTPUT);
		if (this.VERBOSE_BEST_DONOR_OUTPUT != null)
			IOUtil.assertFileIsWritable(this.VERBOSE_BEST_DONOR_OUTPUT);
		// validate VCF is indexed.
		if (!VCFUtils.hasIndex(this.VCF))
			list.add("VCF is not indexed!  Please index and retry.");

		if (this.CELL_CONTAMINATION_ESTIMATE_FILE != null && this.ALLELE_FREQUENCY_ESTIMATE_FILE == null)
			list.add("If CELL_CONTAMINATION_ESTIMATE_FILE is supplied, must also supply ALLELE_FREQUENCY_ESTIMATE_FILE");

		if (this.CELL_CONTAMINATION_ESTIMATE_FILE == null && this.ALLELE_FREQUENCY_ESTIMATE_FILE != null)
			list.add("If ALLELE_FREQUENCY_ESTIMATE_FILE is supplied, must also supply CELL_CONTAMINATION_ESTIMATE_FILE");

		if (this.CELL_CONTAMINATION_ESTIMATE_FILE != null && this.ALLELE_FREQUENCY_ESTIMATE_FILE != null && this.MAX_ERROR_RATE != null) {
			log.info("Running analysis with contamination esimates, MAX_ERROR_RATE argument ignored!");
			MAX_ERROR_RATE = null;
		}

		if (this.CELL_CONTAMINATION_ESTIMATE_FILE != null && this.ALLELE_FREQUENCY_ESTIMATE_FILE != null && this.FIXED_ERROR_RATE != null) {
			log.info("Running analysis with contamination esimates, FIXED_ERROR_RATE argument ignored!");
			FIXED_ERROR_RATE = null;
		}

		if (this.FIXED_ERROR_RATE != null & this.MAX_ERROR_RATE == null)
			log.info("Running with a fixed error rate of " + this.FIXED_ERROR_RATE);

		if (this.MAX_ERROR_RATE != null && this.CELL_CONTAMINATION_ESTIMATE_FILE == null && this.ALLELE_FREQUENCY_ESTIMATE_FILE == null)
			log.info("Running with a maximum cap on error rate of " + this.FIXED_ERROR_RATE);

		if (VERBOSE_BEST_DONOR_OUTPUT != null & VERBOSE_OUTPUT == null) {
			list.add("If VERBOSE_BEST_DONOR_OUTPUT is requested, then VERBOSE_OUTPUT must also be requested");
		}

		if (FRACTION_SAMPLES_PASSING < 0 | FRACTION_SAMPLES_PASSING > 1)
			list.add("FRACTION_SAMPLES_PASSING must be between 0 and 1, value was " + Double.toString(this.FRACTION_SAMPLES_PASSING));
		
		if (DNA_MODE) {
			log.info("IN DNA MODE.  Edit distance for UMI collapse disabled!");
			this.EDIT_DISTANCE=0;
		}

		return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);
	}

}
