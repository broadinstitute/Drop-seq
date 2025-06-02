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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.multisample;

import java.io.File;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.GeneFunctionCommandLineBase;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPUMIBasePileupIterator;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SortOrder;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.CellAssignmentUtils;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.CellCollectionSampleLikelihoodCollection;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.SampleGenotypeProbabilities;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.SampleGenotypeProbabilitiesIterator;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.AssertSequenceDictionaryIntersection;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import org.broadinstitute.dropseqrna.utils.FileUtils;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.VCFUtils;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.IgnoreGeneAnnotationTagger;
import org.broadinstitute.dropseqrna.utils.readiterators.PCRDuplicateFilteringIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.vcftools.SampleAssignmentVCFUtils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import picard.cmdline.StandardOptionDefinitions;
import picard.nio.PicardHtsPath;

@CommandLineProgramProperties(summary = "Detect Doublets in Dropulation Data.  Uses the outputs of AssignCellsToSamples to make decisions.  It's highly recommended to use the VCF output from AssignCellsToSamples as input, as"
		+ "the memory usage of a full VCF may be prohibitive compared to the AssignCellsToSamples VCF, which contains only the variants that were observed in the data.  This also greatly speeds up"
		+ "analysis.", oneLineSummary = "Detect Doublets in Dropulation Data", programGroup = DropSeq.class)
public class DetectDoublets extends GeneFunctionCommandLineBase {

	private static final Log log = Log.getInstance(DetectDoublets.class); 

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze. This argument can accept wildcards, or a file with the suffix .bam_list that contains the locations of multiple BAM files", minElements = 1)
	public List<PicardHtsPath> INPUT_BAM;

	@Argument(doc = "The input VCF file to analyze.  Use the output VCF from AssignCellsToSamples to save memory.")
	public PicardHtsPath VCF;

	@Argument(doc = "The output likelihood file from AssignCellsToSamples")
	public File SINGLE_DONOR_LIKELIHOOD_FILE;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "Output file of doublet likelihoods. This supports zipped formats like gz and bz2.")
	public File OUTPUT;

	@Argument(doc = "Output file of per-sample-pair doublet likelihoods. Insted of just the best pair as seen in the OUTPUT file, "
			+ "this outputs all tested pairs for each cell. This supports zipped formats like gz and bz2.", optional = true)
	public File OUTPUT_ALL_PAIRS = null;

	@Argument(doc = "Output file of per-snp/sample-pair doublet likelihoods. Insted of just the best pair as seen in the OUTPUT file, "
			+ "this outputs all tested pairs for each cell, for each SNP. This file can get pretty obnoxiously huge."
			+ "This supports zipped formats like gz and bz2.", optional = true)
	public File OUTPUT_PER_SNP = null;

	@Argument(doc = "The cell barcode tag.  If there are no reads with this tag, the program will assume that all reads belong to the same cell and process in single sample mode.")
	public String CELL_BARCODE_TAG = "XC";

	@Argument(doc = "The molecular barcode tag.")
	public String MOLECULAR_BARCODE_TAG = "XM";

	@Argument(doc = "The edit distance that molecular barcodes should be combined at within a gene/SNP.")
	public Integer EDIT_DISTANCE = 1;

	@Argument(doc = "The map quality of the read to be included.")
	public Integer READ_MQ = 10;

	@Argument(doc = "Override NUM_CORE_BARCODES and process reads that have the cell barcodes in this file instead.  The file has 1 column with no header.", optional = false)
	public File CELL_BC_FILE = null;

	@Argument(doc = "The minimum genotype quality for a variant.  Set this value to 0 to not filter by GQ scores if they are present, or to -1 to completely "
			+ "ignore GQ values if they are not set in the genotype info field.  If the GQ field is not set in the VCF header, this will be set to -1 by default.")
	public Integer GQ_THRESHOLD = 30;
	
	@Argument(doc = "A file with a list of samples in the VCF to consider as samples in the doublets.  This subsets the VCF into a smaller data set containing only the samples listed. "
			+ "The file has 1 column with no header.  If this list contains only one donor and no contaminating donors were found via single donor assignment, "
			+ "doublet detection calculations will not take place.  Instead, the program will emit a default output for each cell, and the program "
			+ "will then quit with an exit status.", optional = false)
	public File SAMPLE_FILE;

	@Argument(doc = "Instead of useing base qualities to determine error rate, use a fixed error rate instead. This is rounded to the nearest phread score internally.", optional = true)
	public Double FIXED_ERROR_RATE = null;

	@Argument(doc = "Caps the base error rate at a maximum probability so no SNP can be weighed more than this value.  For example, if this value was 0.01, "
			+ "then a base quality 30 value (normally an erro rate of 0.001) would become 0.01.  With the same threshold, a base with an error rate of 0.1 would be unaffected.", optional = true)
	public Double MAX_ERROR_RATE = null;

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

	@Argument(doc = "At least <FRACTION_SAMPLES_PASSING> samples must have genotype scores >= GQ_THRESHOLD for the variant in the VCF to be included in the analysis.")
	public double FRACTION_SAMPLES_PASSING = 0.5;

	@Argument(doc = "A list of chromosomes to omit from the analysis.  The default is to omit the sex chromosomes.")
	public List<String> IGNORED_CHROMOSOMES = new ArrayList<>(Arrays.asList("X", "Y", "MT"));

	/**
	 * this produces worse results with in-silico testing, generating more misclassifications of singlets as doublets.
	 * 
	 * @Argument(doc="For SNPs that don't have a high quality genotype in the VCF, should we infer a global penalty per SNP
	 *                    and apply it to donors for that SNP that are not confidently called?") public boolean
	 *                    USE_MISSING_DATA=true;
	 */
	private boolean USE_MISSING_DATA = false;

	@Argument(doc = "Force evaluation of the doublet at the given mixture ratio.  Should be a number between 0 and 1.", optional = true)
	public Double FORCED_RATIO = 0.8;

	@Argument(doc = "Should cells that were assigned to a donor not on the sample list be tested?  When enabled, "
			+ "this forces the program to load genotype information for all donors seen at least once, even when not on the sample list.  If there are many off-target assignments,"
			+ "this can use large amounts of memory.  Set to false to skip testing cells that you'll probably discard as incorrectly assigned later anyway, set to true to "
			+ "test each cell not on the donor lists against all possible donors on the list.")
	public Boolean TEST_UNEXPECTED_DONORS = true;

	@Argument(doc = "For each cell, when comparing donor pairs to each other, scale the likelihoods to the number of UMIs.  "
			+ "This provides an additional penalty score to donor pairs with more incomplete information that makes donor pairs more comparable.")
	public Boolean SCALE_LIKELIHOODS = true;

	@Argument(doc = "EXPERIMENTAL!!! Run the program in DNA Mode.  In this mode, reads should have a cell barcode, but will be missing gene annotations and UMIs.  All reads will be "
			+ "accepted as passing, and each read (or read pair) will be treated as a single UMI  If the data is PCR Duplicate marked, duplicate reads will be filtered. ")
	public Boolean DNA_MODE = false;

	/*
	 * @Argument
	 * (doc="The model tests the best donor against the all the possible second most likely donors to find the pair that best explain the data.  Sort the donors by their single donor likelihood score,"
	 * +
	 * "and only test the best <MIN_DONOR_PAIRS_TESTED> donors.  This can reduce the number of tests / runtime / memory for large pools while producing approximately the same output."
	 * )
	 */

	private Integer MIN_DONOR_PAIRS_TESTED = Integer.MAX_VALUE;

	// @Argument(doc="See MIN_DONOR_PAIRS_TESTED. This limits the number of donors tested to a fraction of the total number
	// of donors. The number of donors used is the maximum of FRACTION_DONOR_PAIRS_TESTED and MIN_DONOR_PAIRS_TESTED.")

	private Double FRACTION_DONOR_PAIRS_TESTED = 1.0;

	private final String SNP_TAG = "YS";

	private static DecimalFormat mixtureFormat = new DecimalFormat("#.###");

	@Override
	protected int doWork() {
		if (CELL_CONTAMINATION_ESTIMATE_FILE != null) {
			IOUtil.assertFileIsReadable(CELL_CONTAMINATION_ESTIMATE_FILE);
		}
		if (ALLELE_FREQUENCY_ESTIMATE_FILE != null) {
			IOUtil.assertFileIsReadable(ALLELE_FREQUENCY_ESTIMATE_FILE);
		}

		List<String> donorList = ParseBarcodeFile.readCellBarcodeFile(this.SAMPLE_FILE);
		log.info("Number of donors in donor list [" + donorList.size() + "]");
		int pairsToTest = getNumberOfDonorsToTest(donorList.size(), this.MIN_DONOR_PAIRS_TESTED, this.FRACTION_DONOR_PAIRS_TESTED);

		PrintStream perDonorWriter = null;
		if (OUTPUT_ALL_PAIRS != null) {
			perDonorWriter = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(this.OUTPUT_ALL_PAIRS));
			writeAssignmentHeader(perDonorWriter, pairsToTest, false, true);
		}
	
		PrintStream perSNPWriter = null;
		if (OUTPUT_PER_SNP != null) {
			IOUtil.assertFileIsWritable(this.OUTPUT_PER_SNP);
			perSNPWriter = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(this.OUTPUT_PER_SNP));
			writePerSNPHeader(perSNPWriter);
		}

		PrintStream writer = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(this.OUTPUT));
		writeAssignmentHeader(writer, pairsToTest, true, false);

		final VCFFileReader vcfReader = new VCFFileReader(this.VCF.toPath(), false);
		final SamHeaderAndIterator headerAndIter = SamFileMergeUtil.mergeInputPaths(
				PicardHtsPath.toPaths(this.INPUT_BAM), false, SamReaderFactory.makeDefault());
		AssertSequenceDictionaryIntersection.assertIntersectionObjectVcf(
				headerAndIter.header, "BAM INPUT(S)", this.VCF.toPath(), log);

		// extract the sequence dictionary to build the interval list.
		SAMSequenceDictionary vcfDict = vcfReader.getFileHeader().getSequenceDictionary();

		// disable GQ filter if it's not in the header.
		if (!VCFUtils.GQInHeader(vcfReader)) {
			this.GQ_THRESHOLD = -1;
			log.info("Genotype Quality [GQ] not found in header.  Disabling GQ_THRESHOLD parameter");
		}

		CellCollectionSampleLikelihoodCollection cslc = CellCollectionSampleLikelihoodCollection.parseFile(this.SINGLE_DONOR_LIKELIHOOD_FILE);

		// A map where the key is the cell barcode, and the value is the best donor.
		Map<String, String> bestDonorForCell = getBestDonorForCell(this.SINGLE_DONOR_LIKELIHOOD_FILE);

		// Filter the single donor assignment by the donor list if requested.
		if (!this.TEST_UNEXPECTED_DONORS)
			bestDonorForCell = filterDonorMap(bestDonorForCell, donorList);

		// read in the per-cell penalty score and validate against the best donor per cell to make sure you have penalties for
		// every donor.
		Map<String, Double> contaminationMap = CellAssignmentUtils.getCellContamination(this.CELL_CONTAMINATION_ESTIMATE_FILE, bestDonorForCell.keySet());
		Map<Interval, Double> variantMinorAlleleFrequency = CellAssignmentUtils.getMinorAlleleFrequencyMap(this.ALLELE_FREQUENCY_ESTIMATE_FILE);

		// all donors is the input set of donors plus any best donor for a cell.
		// this list of donors is restricted to the sample list if TEST_UNEXPECTED_DONORS=false.
		Set<String> allDonors = new HashSet<>();
		allDonors.addAll(donorList);
		allDonors.addAll(new HashSet<>(bestDonorForCell.values()));
		List<String> allDonorsList = new ArrayList<>(allDonors);

		log.info("Number of donors that can be either donor in a donor pair [" + allDonorsList.size() + "]");

		// Keep all the barcodes that are in the barcode list AND in the single donor assignments.
		List<String> cellBarcodes = getCellBarcodes(this.CELL_BC_FILE, bestDonorForCell, this.TEST_UNEXPECTED_DONORS);
		
		// If there is only one donor at this point, doublet detection should not continue.  
		if (allDonorsList.size()<2) {
			singleDonorGracefulExit(bestDonorForCell, writer, perDonorWriter, perSNPWriter);
			return 0;
		}
		
		// Pass a list of all donors requested + best calls if TEST_UNEXPECTED_DONORS is true.
		// set the %passing to be 0, since the snpIntervals will properly filter on the right set of donors, but you want a
		// super-set of donors available
		// so if you call donors A-D, but the single cell assigned E, then E would still be an available donor in the genotype
		// matrix.
		PeekableIterator<VariantContext> vcfIterator = SampleAssignmentVCFUtils.getVCFIterator(vcfReader, allDonorsList, false, this.GQ_THRESHOLD,
				this.FRACTION_SAMPLES_PASSING, IGNORED_CHROMOSOMES, log);
		
		GenotypeMatrix genotypeMatrix = new GenotypeMatrix(vcfIterator, this.GQ_THRESHOLD, allDonorsList);
		vcfIterator.close();

		Map<Interval, Double> genotypeQuality = genotypeMatrix.getAverageGenotypeQuality();
		
		final IntervalList snpIntervals = new IntervalList(vcfDict);								
		snpIntervals.addall(genotypeMatrix.getSNPIntervals());
		// a requested early exit if there are no SNPs.
		if (snpIntervals.getIntervals().isEmpty()) {
			log.error("No SNP intervals detected!  Check to see if your VCF filter thresholds are too restrictive!");
			return 1;
		}
				
		PeekableIterator<List<SampleGenotypeProbabilities>> sampleGenotypeIterator = prepareIterator(snpIntervals, cellBarcodes, genotypeQuality);

		int cellCount = 0;
		int reportInterval = 100;
		log.info("Calling doublets");
		if (!sampleGenotypeIterator.hasNext()) {
			log.warn("No Cells found for analysis.");
		} else {

			while (sampleGenotypeIterator.hasNext()) {
				cellCount++;
				if (cellCount % reportInterval == 0)
					log.info("Tested cell #" + cellCount);
				List<SampleGenotypeProbabilities> probs = sampleGenotypeIterator.next();

				String cell = probs.get(0).getCell();
				String bestDonor = bestDonorForCell.get(cell);
				if (bestDonor == null)
					throw new IllegalStateException("Cell [" + cell + "] has no best donor assignment in file.");

				VariantDataFactory f = null;
				f = new VariantDataFactory(cell, probs, genotypeMatrix, FIXED_ERROR_RATE, USE_MISSING_DATA, MAX_ERROR_RATE, contaminationMap,
						variantMinorAlleleFrequency);

				FindOptimalDonorMixture fodm = new FindOptimalDonorMixture(f);

				// AllPairedSampleAssignmentsForCell allAssignments = fodm.findBestDonorPair(bestDonor, donorList, FORCED_RATIO);
				// List<String> donorsThisCell = getExpectedSecondDonorsRankedByLikelihood(cell, cslc, pairsToTest, donorList, bestDonor);
				AllPairedSampleAssignmentsForCell allAssignments = fodm.findBestDonorPair(bestDonor, donorList, FORCED_RATIO, SCALE_LIKELIHOODS);
				SamplePairAssignmentForCell best = allAssignments.getBestAssignment();

				// edge case: assignment is null because there's no data for this cell.  This only happens
				// when a user has a cell selection error or similar and attempts to call cell barcodes that aren't cells.				
				if (best==null) { 
					log.warn("No best pair found for cell ["+cell+"] due to no informative UMIs.  Cell selection or similar problem?");
					best= SamplePairAssignmentForCell.constructEmptyResult(cell, bestDonorForCell.get(cell));
				}
				 
				double bestPairPvalue = allAssignments.getBestPairPvalue();																	
				writeAssignment(best, bestPairPvalue, writer, false);
				if (OUTPUT_ALL_PAIRS != null) {
					writeAssignment(allAssignments.getBestAssignment(), null, perDonorWriter, true);
					// apply ordering to other assignments for stability of outputs. Sort by 2nd donor name.
					Comparator<SamplePairAssignmentForCell> comparator = java.util.Comparator.comparing(SamplePairAssignmentForCell::getSampleTwo,
							java.util.Comparator.naturalOrder());
					List<SamplePairAssignmentForCell> allOther = allAssignments.getOtherAssignments();
					Collections.sort(allOther, comparator);
					for (SamplePairAssignmentForCell other : allOther)
						writeAssignment(other, null, perDonorWriter, true);
				}
				reportResultsPerSNP(cell, f, bestDonor, donorList, allAssignments, perSNPWriter);
			}
		}
		
		if (OUTPUT_PER_SNP != null)
			perSNPWriter.close();
		if (OUTPUT_ALL_PAIRS != null)
			perDonorWriter.close();
		writer.close();

		log.info("Finished!");
		return 0;
	}
	
	
	/**
	 * In the strange edge case where there is only a single donor to be tested, write out a default output file instead of going through testing.
	 * This function additionally closes all potentially open writers and runs logging.
	 * @param bestDonorForCell A map containing cell barcodes and the best donor for each cell.
	 * @param writer The file to write to per-cell outputs.	 
	 * @param perDonorWriter Closes this file if not null.
	 * @param perSNPWriter Closes this file if not null.
	 */
	void singleDonorGracefulExit(Map<String, String> bestDonorForCell, PrintStream writer, 
			PrintStream perDonorWriter, PrintStream perSNPWriter) {
		// clean up more detailed file writers. 
		if (OUTPUT_ALL_PAIRS!=null) perDonorWriter.close();
		if (OUTPUT_PER_SNP != null) perSNPWriter.close();
		// write a default output per cell close results and quit.			
		log.error("The donor file only contained a single donor, and no additional donors were detected by single donor assignment.  Doublet detection will not continue.  "
				+ "A default output will be written to perserve downstream pipeline functionality.");
		writeSingleDonorEdgeCaseOutput(bestDonorForCell, writer);				
	}
	
	/**
	 * In the strange edge case where there is only a single donor to be tested, write out a default output file instead of going through testing.
	 * @param bestDonorForCell A map containing cell barcodes and the best donor for each cell.
	 * @param writer The file to write to
	 */
	void writeSingleDonorEdgeCaseOutput(Map<String, String> bestDonorForCell, PrintStream writer) {
		for (String cell: bestDonorForCell.keySet()) {
			SamplePairAssignmentForCell best = SamplePairAssignmentForCell.constructEmptyResult(cell, bestDonorForCell.get(cell));
			writeAssignment(best, 1.0E-101d, writer, false);
		}
		writer.close();
		
	}

	/**
	 * Filter the map of cell barcode -> donor to only retain cell barcodes where the assigned donor is in the donor list.
	 * 
	 * @param bestDonorForCell
	 * @param donorList
	 * @return A submap of the input map where all values are contained in the donor list.
	 */
	Map<String, String> filterDonorMap(Map<String, String> bestDonorForCell, List<String> donorList) {
		Set<String> dl = new HashSet<String>(donorList);
		Map<String, String> result = new HashMap<String, String>();

		for (String cellBC : bestDonorForCell.keySet()) {
			String donor = bestDonorForCell.get(cellBC);
			if (dl.contains(donor))
				result.put(cellBC, donor);
		}
		return result;
	}

	private void reportResultsPerSNP(final String cell, final VariantDataFactory variantFactory, final String bestDonor, final List<String> vcfSamples,
			final AllPairedSampleAssignmentsForCell allAssignments, final PrintStream out) {
		if (out == null)
			return;
		List<String> other = FindOptimalDonorMixture.getNonPrimarySamples(bestDonor, vcfSamples);
		for (String o : other) {
			VariantDataCollection vdc = variantFactory.getVariantData(bestDonor, o);
			SamplePairAssignmentForCell mixtureResult = allAssignments.getAssignmentForDonorPair(bestDonor, o);
			List<VariantData> vdList = vdc.getVariantData();
			double mixture = mixtureResult.getMixture();
			for (VariantData vd : vdList)
				writePerSNPReport(cell, vd, bestDonor, o, mixture, out);
		}
	}

	private void writePerSNPReport(final String cell, final VariantData vd, final String sampleOne, final String sampleTwo, final Double mixture,
			final PrintStream out) {
		/*
		 * if (mixture==null) { String [] line = {cell, sampleOne, sampleTwo, vd.getSNPInterval().getContig(),
		 * Integer.toString(vd.getSNPInterval().getStart()), vd.getGenotypeOne().toString(), vd.getGenotypeTwo().toString(),
		 * Integer.toString(vd.getGenotypeCountReference()), Integer.toString(vd.getGenotypeCountAlternate()), "NA",
		 * Double.toString(vd.getLogLikelihood(1)), Double.toString(vd.getLogLikelihood(0))}; String h = StringUtils.join(line,
		 * "\t"); out.println(h); return; }
		 */

		String[] line = { cell, sampleOne, sampleTwo, vd.getSNPInterval().getContig(), Integer.toString(vd.getSNPInterval().getStart()),
				vd.getGenotypeOne().toString(), vd.getGenotypeTwo().toString(), Integer.toString(vd.getGenotypeCountReference()),
				Integer.toString(vd.getGenotypeCountAlternate()), Double.toString(vd.getLogLikelihood(mixture)), Double.toString(vd.getLogLikelihood(1)),
				Double.toString(vd.getLogLikelihood(0)) };
		String h = StringUtils.join(line, "\t");
		out.println(h);
	}

	private void writePerSNPHeader(final PrintStream out) {
		String[] line = { "cell", "sampleOne", "sampleTwo", "chr", "pos", "genotype_S1", "genotype_S2", "refAlleleCount", "altAlleleCount",
				"likelihood_mixture", "likelihood_S1", "likelihood_S2" };
		String h = StringUtils.join(line, "\t");
		out.println(h);
	}

	private String convertNullToString(final Double x) {
		if (x == null)
			return ("NA");
		return Double.toString(x);
	}

	private void writeAssignmentHeader(final PrintStream out, final int pairsToTest, final boolean outputBestPairPvalue, final boolean writeScaledLikelihoods) {
		final List<String> paths = FileUtils.toAbsoluteStrings(PicardHtsPath.toPaths(this.INPUT_BAM));
		String bamList = StringUtils.join(paths, ",");

		final List<String> header = new ArrayList<>(Arrays.asList(
				"#INPUT_BAM=" + bamList, "INPUT_VCF=" + FileUtils.toAbsoluteString(this.VCF.toPath()),
				"DONOR_FILE=" + this.SAMPLE_FILE, "CELL_BC_FILE=" + CELL_BC_FILE, "GQ_THRESHOLD=" + Integer.toString(this.GQ_THRESHOLD),
				"FRACTION_SAMPLES_PASSING=" + Double.toString(this.FRACTION_SAMPLES_PASSING), "FORCED_RATIO=" + convertNullToString(FORCED_RATIO),
				"USE_MISSING_DATA=" + USE_MISSING_DATA, "READ_MQ=" + Integer.toString(this.READ_MQ),
				"FIXED_ERROR_RATE=" + convertNullToString(this.FIXED_ERROR_RATE), "MAX_ERROR_RATE=" + convertNullToString(this.MAX_ERROR_RATE),
				"LOCUS_FUNCTION=" + this.LOCUS_FUNCTION_LIST.toString(), "PAIRS_TO_TEST=" + Integer.toString(pairsToTest)));

		if (this.CELL_CONTAMINATION_ESTIMATE_FILE != null && this.ALLELE_FREQUENCY_ESTIMATE_FILE != null) {
			header.add("CELL_CONTAMINATION_ESTIMATE_FILE=" + CELL_CONTAMINATION_ESTIMATE_FILE.getAbsolutePath());
			header.add("ALLELE_FREQUENCY_ESTIMATE_FILE=" + ALLELE_FREQUENCY_ESTIMATE_FILE.getAbsolutePath());
		}

		String h = StringUtils.join(header, "\t");
		out.println(h);

		writeAssignmentColumnNames(out, outputBestPairPvalue, writeScaledLikelihoods);
	}

	public static void writeAssignmentColumnNames(final PrintStream out, final boolean outputBestPairPvalue, final boolean writeScaledLikelihoods) {
		List<String> line = new ArrayList<String>(Arrays.asList("cell", "sampleOneMixtureRatio", "sampleOne", "sampleOneLikelihood", "sampleTwo",
				"sampleTwoLikelihood", "mixedSample", "mixedSampleLikelihood", "num_paired_snps", "num_inform_snps", "num_umi", "num_inform_umis",
				"lr_test_stat", "sampleOneWrongAlleleCount", "num_homozygous_inform_umis_s1", 
				"sampleTwoWrongAlleleCount", "num_homozygous_inform_umis_s2", "bestLikelihood", "bestSample", "doublet_pval"));

		if (outputBestPairPvalue)
			line.add("best_pair_pvalue");

		if (writeScaledLikelihoods) {
			line.add("bestLikelihoodScaled");
		}

		String header = StringUtils.join(line, "\t");
		out.println(header);
	}

	public static void writeAssignment(final SamplePairAssignmentForCell assignment, final Double bestPairPvalue, final PrintStream out,
			final boolean writeScaledLikelihood) {				
		String mixture = mixtureFormat.format(assignment.getMixture());

		List<String> line = new ArrayList<String>(
				Arrays.asList(assignment.getCellBarcode(), mixture, assignment.getSampleOne(), Double.toString(assignment.getSampleOneSingleLikelihood()),
						assignment.getSampleTwo(), Double.toString(assignment.getSampleTwoSingleLikelihood()), assignment.getCombinedDonorName(),
						Double.toString(assignment.getDoubletLikelihood()), Integer.toString(assignment.getNumSNPs()),
						Integer.toString(assignment.getNumInformativeSNPs()), Integer.toString(assignment.getNumUMIs()),
						Integer.toString(assignment.getNumInformativeUMIs()), Double.toString(assignment.getDoubletLikelihoodRatio()),
						Integer.toString(assignment.getImpossibleAllelesSampleOne()), Integer.toString(assignment.getNumInformativeHomozygousUMIsSampleOne()), 
						Integer.toString(assignment.getImpossibleAllelesSampleTwo()), Integer.toString(assignment.getNumInformativeHomozygousUMIsSampleTwo()),
						Double.toString(assignment.getBestLikelihood()), assignment.getBestSample(), Double.toString(assignment.getDoubletPvalue())));

		if (bestPairPvalue != null)
			line.add(bestPairPvalue.toString());
		
		if (writeScaledLikelihood)
			line.add(Double.toString(assignment.getScaledBestLikelihood()));

		String h = StringUtils.join(line, "\t");
		out.println(h);
	}
	
	public PeekableIterator<List<SampleGenotypeProbabilities>> prepareIterator(final IntervalList snpIntervals, List<String> cellBarcodes, Map<Interval, Double> genotypeQuality) {

		SamReaderFactory factory = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE);
		SamHeaderAndIterator headerAndIter =
				SamFileMergeUtil.mergeInputPaths(PicardHtsPath.toPaths(this.INPUT_BAM), false, factory);

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
		
		SNPUMIBasePileupIterator sbpi = new SNPUMIBasePileupIterator(headerAndIter, snpIntervals, GENE_NAME_TAG, GENE_STRAND_TAG, GENE_FUNCTION_TAG,
				LOCUS_FUNCTION_LIST, STRAND_STRATEGY, this.FUNCTIONAL_STRATEGY, this.CELL_BARCODE_TAG, this.MOLECULAR_BARCODE_TAG, this.SNP_TAG,
				GeneFunctionCommandLineBase.DEFAULT_FUNCTION_TAG, this.READ_MQ, false, cellBarcodes, genotypeQuality, SortOrder.CELL_SNP);

		final SAMSequenceDictionary dict = snpIntervals.getHeader().getSequenceDictionary();

		// gets a SampleGenotypeProbabilities for each cell.
		SampleGenotypeProbabilitiesIterator result = new SampleGenotypeProbabilitiesIterator(sbpi, dict, this.EDIT_DISTANCE, SortOrder.CELL_SNP);

		// clusters SampleGenotypeProbabilities objects across all cells for a SNP
		GroupingIterator<SampleGenotypeProbabilities> groupingIterator = new GroupingIterator<>(result, new Comparator<SampleGenotypeProbabilities>() {
			@Override
			public int compare(final SampleGenotypeProbabilities o1, final SampleGenotypeProbabilities o2) {
				int cmp = o1.getCell().compareTo(o2.getCell());
				return cmp;
			}
		});

		PeekableIterator<List<SampleGenotypeProbabilities>> peekableIter = new PeekableIterator<>(groupingIterator);
		return (peekableIter);
	}

	List<String> getCellBarcodes(File cellBarcodeFile, Map<String, String> bestDonorForCell, boolean testUnexpectedDonors) {
		List<String> cellBarcodes = ParseBarcodeFile.readCellBarcodeFile(this.CELL_BC_FILE);
		int numBarcodes = cellBarcodes.size();

		log.info("Number of cell barcodes in input file [" + numBarcodes + "]");

		if (!testUnexpectedDonors) {
			cellBarcodes.retainAll(bestDonorForCell.keySet());
			log.info("Number of cell barcodes after filtering to expected donors [" + cellBarcodes.size() + "]");
			double fracRemoved = (double) (numBarcodes - cellBarcodes.size()) / (double) numBarcodes;
			log.info("% cell barcodes not assigned to an expected donor [" + new DecimalFormat("0.##").format(fracRemoved * 100) + "%]");
		}

		return cellBarcodes;
	}

	/**
	 * Parse the single donor likelihood file, retrieve the best donor for each cell barcode.
	 * 
	 * @param singleDonorLikelihoodFile The OUTPUT file produced by AssignCellsToSamples.
	 * @return A map where the key is the cell barcode, and the value is the best donor.
	 */
	Map<String, String> getBestDonorForCell(final File singleDonorLikelihoodFile) {
		Map<String, String> result = new HashMap<>();

		CellCollectionSampleLikelihoodCollection cslc = CellCollectionSampleLikelihoodCollection.parseFile(singleDonorLikelihoodFile);

		for (String cellBarcode : cslc.getCellBarcodes()) {
			String donor = cslc.getLikelihoodCollection(cellBarcode).getBestSampleAssignment().getSample();
			result.put(cellBarcode, donor);
		}
		return (result);
	}

	// getNumberOfDonorsToTest(donorList.size(), this.MIN_DONOR_PAIRS_TESTED, this.FRACTION_DONOR_PAIRS_TESTED);
	private int getNumberOfDonorsToTest(int totalDonorPairs, int minDonorPairs, double fractionDonorPairs) {
		int numTest = (int) Math.ceil((double) totalDonorPairs * fractionDonorPairs);
		int result = Math.max(numTest, minDonorPairs);
		// don't test more than the total number of donors.
		// this is for when total donors are < minDonorPairs.
		if (result > totalDonorPairs)
			result = totalDonorPairs;
		log.info("Testing " + Integer.toString(result) + " Donor pairs per cell");
		return (result);

	}

		
	@Override
	protected String[] customCommandLineValidation() {
		IOUtil.assertFileIsReadable(this.VCF.toPath());
		this.INPUT_BAM = FileListParsingUtils.expandPicardHtsPathList(INPUT_BAM);
		IOUtil.assertFileIsWritable(this.OUTPUT);
		IOUtil.assertFileIsReadable(this.SAMPLE_FILE);
		IOUtil.assertFileIsReadable(this.SINGLE_DONOR_LIKELIHOOD_FILE);
		IOUtil.assertFileIsReadable(this.CELL_BC_FILE);
		final ArrayList<String> list = new ArrayList<>(1);
		if (OUTPUT_ALL_PAIRS != null)
			IOUtil.assertFileIsWritable(this.OUTPUT_ALL_PAIRS);
		if (OUTPUT_PER_SNP != null)
			IOUtil.assertFileIsWritable(this.OUTPUT_PER_SNP);

		if (!VCFUtils.hasIndex(this.VCF.toPath()))
			list.add("VCF is not indexed!  Please index and retry.");

		if (this.CELL_CONTAMINATION_ESTIMATE_FILE != null && this.ALLELE_FREQUENCY_ESTIMATE_FILE == null)
			list.add("If CELL_CONTAMINATION_ESTIMATE_FILE is supplied, must also supply ALLELE_FREQUENCY_ESTIMATE_FILE");

		if (this.CELL_CONTAMINATION_ESTIMATE_FILE == null && this.ALLELE_FREQUENCY_ESTIMATE_FILE != null)
			list.add("If ALLELE_FREQUENCY_ESTIMATE_FILE is supplied, must also supply CELL_CONTAMINATION_ESTIMATE_FILE");

		if (this.FIXED_ERROR_RATE != null & this.MAX_ERROR_RATE == null)
			log.info("Running with a fixed error rate of " + this.FIXED_ERROR_RATE);

		if (this.MAX_ERROR_RATE != null && this.CELL_CONTAMINATION_ESTIMATE_FILE == null && this.ALLELE_FREQUENCY_ESTIMATE_FILE == null)
			log.info("Running with a maximum cap on error rate of " + this.MAX_ERROR_RATE);

		if (FRACTION_SAMPLES_PASSING < 0 | FRACTION_SAMPLES_PASSING > 1)
			list.add("FRACTION_SAMPLES_PASSING must be between 0 and 1, value was " + Double.toString(this.FRACTION_SAMPLES_PASSING));

		return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), list);
	}

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new DetectDoublets().instanceMain(args));
	}
	
	/*
	private List<String> getExpectedSecondDonorsRankedByLikelihood(String cellBarcode, CellCollectionSampleLikelihoodCollection cslc, int pairsToTest,
			List<String> expectedDonors, String bestDonor) {
		// short circuit if the requested number of donors is the same as the number of expected donors IE the non-optimized
		// strategy.
		if (expectedDonors.size() == pairsToTest)
			return expectedDonors;
		
		CellSampleLikelihoodCollection c = cslc.getLikelihoodCollection(cellBarcode);
		if (c == null)
			throw new IllegalArgumentException("Could not find single donor likelihoods for cell " + cellBarcode);
		List<String> rankedDonors = c.getDonorsRankedByAssignmentLikelihood();
		// exclude the best donor from the ranked list. We don't want to select that.
		rankedDonors.remove(bestDonor);

		Set<String> expected = new HashSet<String>(expectedDonors);
		// filter ranked donors by expected.
		List<String> rankedExpectedDonors = rankedDonors.stream().filter(x -> expected.contains(x)).collect(Collectors.toList());

		// if pairs to test is set too high somehow (the sample lists doesn't match up with the input single donor assignments)
		// then limit the return.
		if (rankedExpectedDonors.size() < pairsToTest) {
			pairsToTest = rankedExpectedDonors.size();
		}
		// get the top <X> donors.
		rankedExpectedDonors = rankedExpectedDonors.subList(0, pairsToTest);
		return rankedExpectedDonors;

	}
	*/


}
