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
package org.broadinstitute.dropseqrna.censusseq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.censusseq.VCFPileupJointIterator.JointResult;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.AssertSequenceDictionaryIntersection;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.VCFUtils;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.utils.statistics.Diversity;
import org.broadinstitute.dropseqrna.vcftools.SampleAssignmentVCFUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(summary = "Given a VCF file, a list of donors from that file, and BAM file, find the relative fraction of each donor contained in the sequencing experiment.  Handles both single and paired end reads. By default this filters reads with the PCR Duplicate flag.", oneLineSummary = "Create singleton SNP VCF", programGroup = DropSeq.class)
public class CensusSeq extends CommandLineProgram {

	private static final Log log = Log.getInstance(CensusSeq.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze. This argument can accept wildcards, or a file with the suffix .bam_list that contains the locations of multiple BAM files", minElements = 1)
	public List<File> INPUT_BAM;
	
	@Argument(doc = "The input VCF.")
	public File INPUT_VCF;

	@Argument (doc="A file containing a list of samples to filter the VCF by.  Has 1 column, 1 entry per row.  Each entry is a single sample name from the VCF.  If this list is not specified, then all donors in the VCF are used.", optional=true)
	public File SAMPLE_FILE;

	@Argument (shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output file containing each donor in the input SAMPLE_FILE and the relative fraction of the donor.  When run in private SNP mode, an additional column REP_IRVs is emitted.  This is the fraction of the identifying rare alleles - alleles that private to the donor and should only be observed when the donor is present.", optional=false)
	public File OUTPUT;

	@Argument (doc="Output a version of the VCF containing only the variants that have coverage in the input BAM.  This file can be used as the VCF input file on subsequent runs of the same data set to greatly speed up the run, if the same BAM is being analyzed with different conditions.", optional=true)
	public File VCF_OUTPUT;

	@Argument (doc="Output a coverage report for the number of reads on each SNP.  Bins SNPs by the number of reads.", optional=true)
	public File SNP_COVERAGE_HISTOGRAM;

	@Argument(doc= "The minimum genotype quality for a variant.")
	public Integer GQ_THRESHOLD=30;

	@Argument (doc="At least <FRACTION_SAMPLES_PASSING> samples must have genotype scores >= GQ_THRESHOLD for the variant in the VCF to be included in the analysis.")
	public double FRACTION_SAMPLES_PASSING=0.9;

	@Argument (doc="At least this many samples must have the non-common genotype.")
	public Integer MIN_NUM_VARIANT_SAMPLES=2;

	@SuppressWarnings({ "unchecked", "rawtypes" })
	@Argument (doc="A list of chromosomes to omit from the analysis.  The default is to omit the sex chromosomes.")
	public List<String> IGNORED_CHROMOSOMES= new ArrayList(Arrays.asList("X", "Y", "MT"));

	@Argument(doc = "The map quality of the read to be included.")
	public Integer READ_MQ = 10;
	
	@Argument(doc = "The minimum base quality of the SNP base on the read", optional = true)
	public Integer MIN_BASE_QUALITY = 10;

	@Argument(doc = "This option is intended to make evaluation of the algorithm easier when supplied proper training data.  "
			+ "When tagging reads from known donors with the KNOWN_DONOR_TAG and synthetically mixing them to evalutate the algorithm, "
			+ "this counts the number of reads observed with the tag per donor.  Non-informative reads (low map quality, PCR duplicates) are not counted. "
			+ "Reads are only counted on contigs for which there is at least one SNP present in the VCF (minus IGNORED_CHROMOSOES that are ignored) "
			+ "that could inform the data. This tends to exclude bait contigs and alternative haplotypes.", optional=true)
	public String KNOWN_DONOR_TAG=null;

	@Argument (doc="Report how many ref/het/alt/missing genotypes each donor has.", optional=true)
	public Boolean REPORT_ALLELE_COUNTS=false;

	private final String SNP_TAG = "YS";

	@Argument(doc="The number of threads to use for likelihood calculations.")
	public Integer NUM_THREADS=1;

    @Argument(doc = "Random seed to use if deterministic behavior is desired.  " +
            "Setting to null will cause multiple invocations to produce different results.")
    public Integer RANDOM_SEED = 1;

	// How many times will the likelihood try to converge?
	private Integer MAX_RESTARTS=10;

	@Argument(doc="If a SNP has an unusual number of reads, skip it.", optional=true)
	public Integer MAX_READS_PER_SNP=null;

	// @Argument(doc="If convergence doesn't occur, use randomized donor frequency start.")
	// public Boolean RANDOMIZED_START=false;

	@Override
	public int doWork() {
		
		this.INPUT_BAM = FileListParsingUtils.expandFileList(INPUT_BAM);
		IOUtil.assertFileIsReadable(INPUT_VCF);

		// validate VCF is indexed.
		if (!VCFUtils.hasIndex(this.INPUT_VCF)) return 1;

		final VCFFileReader vcfReader = new VCFFileReader(this.INPUT_VCF, true);
		
		SamHeaderAndIterator headerAndIter= SamFileMergeUtil.mergeInputs(this.INPUT_BAM, false, SamReaderFactory.makeDefault());				
		AssertSequenceDictionaryIntersection.assertIntersectionObjectVcf(headerAndIter.header, "BAM INPUT(S)", this.INPUT_VCF, log);		
		SAMSequenceDictionary sd = vcfReader.getFileHeader().getSequenceDictionary();

		if (!VCFUtils.GQInHeader(vcfReader)) {
			this.GQ_THRESHOLD=-1;
			log.info("Genotype Quality [GQ] not found in header.  Disabling GQ_THRESHOLD parameter");
		}
				
		// get the list of samples.  Start with the list of all samples, filter to the sample list if needed.  
		List<String> vcfSamples  = CensusSeqUtils.getFinalSamplelist (vcfReader, this.SAMPLE_FILE);
		if (vcfSamples==null) return 1; // exit if the vcf sample list is null
		
		// write VCF to a temp directory on the first pass so we can iterate on it along with the sorted BAM.
		File outTempVCF=CensusSeqUtils.getTempVCFFile (this.TMP_DIR.get(0));
		VariantContextWriter vcfWriter= CensusSeqUtils.getVCFWriter (outTempVCF, vcfReader, vcfSamples);
					
		// open iterator and scan VCF once to get a list of sites.
		log.info("Looking through VCF for SNPs that fit criteria.  Will search for these in BAM.");
		
		// this filters to SNPs in the sample set.		 
		PeekableIterator<VariantContext> vcfIterator = CensusSeqUtils.getVCFIterator (this.INPUT_VCF, vcfSamples, this.MIN_NUM_VARIANT_SAMPLES, this.GQ_THRESHOLD, 
				this.FRACTION_SAMPLES_PASSING, this.IGNORED_CHROMOSOMES, log, false);
		
		// last argument is a bit funky.  If you aren't using the common SNP analysis, you need to preserve the full interval names.
		final IntervalList snpIntervals = SampleAssignmentVCFUtils.getSNPIntervals(vcfIterator, sd, log, vcfWriter, false);

		if (snpIntervals.size()==0) {
			log.error("0 SNPs found.  Something is very wrong!");
			return 1;
		}
		vcfIterator.close();
		vcfWriter.close();

		// set up the summary output
		BufferedWriter outSummary = null;
		if (this.OUTPUT!=null) {
			IOUtil.assertFileIsWritable(this.OUTPUT);
			outSummary=new BufferedWriter(new OutputStreamWriter(IOUtil.openFileForWriting(OUTPUT)));
			writeDonorRatioOutputHeader(outSummary);
		}

		// build the BAM iterator.
		log.info("Finding SNPs in BAM.");		
		SNPGenomicBasePileupIterator pileUpIter = new SNPGenomicBasePileupIterator(headerAndIter, snpIntervals, SNP_TAG, READ_MQ, this.IGNORED_CHROMOSOMES, KNOWN_DONOR_TAG, this.MIN_BASE_QUALITY);

		// reset the iterator for use in the full data set.  Use the cleaned up set of variants, which should be smaller and faster.
		vcfIterator = CensusSeqUtils.getVCFIterator (outTempVCF, vcfSamples, this.MIN_NUM_VARIANT_SAMPLES, this.GQ_THRESHOLD, 
				this.FRACTION_SAMPLES_PASSING, this.IGNORED_CHROMOSOMES, log, false);
				
		// set up the VCF writer for the final pass, if the VCF_OUTPUT is not null.
		vcfWriter = CensusSeqUtils.getVCFWriter (this.VCF_OUTPUT, vcfReader, vcfSamples);
		
		int runResult=runCensusSeq(vcfSamples, vcfIterator, pileUpIter, snpIntervals.size(), vcfWriter, sd, outSummary, SNP_COVERAGE_HISTOGRAM);		

		vcfReader.close();
		return runResult;
	}


	private int runCensusSeq(final List<String> donorNames, final PeekableIterator<VariantContext> vcfIterator, final SNPGenomicBasePileupIterator pileUpIter,
			final int numPossibleSNPs, final VariantContextWriter vcfWriter, final SAMSequenceDictionary sd, final BufferedWriter outSummary, final File snpHistogramFile) {

		PeekableIterator<SNPGenomicBasePileUp> peekablePileUpIter = new PeekableIterator<>(pileUpIter);

		CommonSNPsData data = new CommonSNPsData(donorNames);

		VCFPileupJointIterator jointIter = new VCFPileupJointIterator(peekablePileUpIter, vcfIterator, sd);
		int numSNPsSkipped=0;
		ObjectCounter<Integer> snpReadDepth = new ObjectCounter<>();
		while (jointIter.hasNext()) {
			JointResult jr = jointIter.next();
			int snpNumBases = jr.getPileup().getNumBases();
			// a bit of a hack.  If the SNP pileup has too many reads, discard it. If parameter not set keep it no matter what.
			if (this.MAX_READS_PER_SNP==null || snpNumBases<=this.MAX_READS_PER_SNP) {
				data = processPileUp (data, jr.getVc(), jr.getPileup(), vcfWriter);
				snpReadDepth.increment(snpNumBases);
			}
			else
				numSNPsSkipped++;
		}

		if (this.MAX_READS_PER_SNP!=null)
			log.info("Number of SNPs skipped [" + numSNPsSkipped+ "] with read count higher than [" + this.MAX_READS_PER_SNP +"]");

		jointIter.close();
		JointIteratorCounter counter = jointIter.getCounter();

		if (counter.SAMPLE_GENE_ITER!=counter.BOTH)
			log.info("Expected counts don't match: Pileup [" + counter.SAMPLE_GENE_ITER + "] both counter [" + counter.BOTH+ "]");

		if (vcfWriter!=null) vcfWriter.close();
		log.info("Optimizing donor ratios across " + data.getNumVariants() +" SNPs and " + donorNames.size() +" donors");
		if (data.getNumVariants()==0 || data.getSampleNames().size()==0) {
			log.error("Number of variants or number of donors is 0.  Can't optimize!");
			return 1;
		}

		// initialize random seed
		final Random random = this.RANDOM_SEED == null ? new Random() : new Random(this.RANDOM_SEED);

		OptimizeSampleRatiosCommonSNPs optim = new OptimizeSampleRatiosCommonSNPs(data, NUM_THREADS, random, false);
		OptimizeSampleRatiosCommonSNPsResult result = optim.directIteration();
		int restartCounter=0;
		while (!result.isConverged() && restartCounter<this.MAX_RESTARTS) {
			restartCounter++;
			log.info("Likelihoods did not converge, restarting with randomized start positions, retry ["+ restartCounter+"].");
			optim = new OptimizeSampleRatiosCommonSNPs(data, NUM_THREADS, random, true);
			OptimizeSampleRatiosCommonSNPsResult resultNew = optim.directIteration();
			if (resultNew.getBestLikelihood()> result.getBestLikelihood()) {
				log.info("Restart with random positions found better maximum likelihood [" + resultNew.getBestLikelihood() +"] old result [" + result.getBestLikelihood()+"]");
				result=resultNew;
			} // if the likelihood converged, you're done.
			else log.info("Random restart likelihood ["+ resultNew.getBestLikelihood()+"] worse than a previous run ["+result.getBestLikelihood()+"].  Rejecting restart result.");
		}

		ObjectCounter<String> knownDonorCounts = pileUpIter.getKnownDonorCounts();

		Map<String, Double> donorRatios = result.getResult();
		writeCommonDonorRatioOutputAdditionalData (outSummary, numPossibleSNPs, counter.BOTH, result, pileUpIter.getKnownDonorTag());
		writeCommonDonorRatioOutput (donorRatios, knownDonorCounts, data, outSummary);
		log.info("Processed [" + (counter.BOTH - counter.REJECTED) + "] SNPs in BAM + VCF");

		writeSNPCountHistogram(snpReadDepth, this.MAX_READS_PER_SNP, numSNPsSkipped, snpHistogramFile);

		if (!result.isConverged()) {
			log.error("Donor ratios did not converge!");
			// if data did not converge, return NOT OK.
			return 1;
		}
		return 0;
	}

	private void writeSNPCountHistogram (final ObjectCounter<Integer> snpCounts, final Integer maxReadsPerSNP, final int numSNPsSkipped, final File snpHistogramFile) {
		List<Integer> keys = new ArrayList<>(snpCounts.getKeys());
		Collections.sort(keys);

		List<String> headerFlex = new ArrayList<>();

		headerFlex.add("#MAX_READS_PER_SNP="+maxReadsPerSNP);
		headerFlex.add("SKIPPED_SNPS_READ_DEPTH="+numSNPsSkipped);

		String [] header = {"READ_DEPTH", "COUNT"};
		if (snpHistogramFile==null) {
			log.info("SNP Read Depth Histogram");
			log.info(StringUtil.join("\t", header));
			for (Integer k: keys)
				log.info(k+"\t"+snpCounts.getCountForKey(k));
		} else {
			IOUtil.assertFileIsWritable(snpHistogramFile);
			PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(snpHistogramFile));
			out.println(StringUtil.join("\t", headerFlex));
			out.println(StringUtil.join("\t", header));
			for (Integer k: keys) {
				String [] line = {Integer.toString(k),Integer.toString(snpCounts.getCountForKey(k))};
				out.println(StringUtil.join("\t", line));
			}
			CloserUtil.close(out);
		}
	}

	private void writeCommonDonorRatioOutputAdditionalData (final BufferedWriter outTerse, final int numPossibleSNPs, final int numSNPsFound, final OptimizeSampleRatiosCommonSNPsResult result, final String knownDonorTag) {
		DecimalFormat divFormat = new DecimalFormat("0.##");

		double [] ratios = result.getResult().values().stream().mapToDouble(x -> x).toArray();
		double diversity = Diversity.diversity(ratios);
		double equitability = Diversity.equitability(ratios);
		String kdt= knownDonorTag;
		if (kdt==null) kdt= "NULL";

		String [] header={"NUM_POSSIBLE_SNPS="+numPossibleSNPs, "NUM_SNPS_USED="+numSNPsFound, "CONVERGED="+result.isConverged(), "BEST_LIKELIHOOD="+result.getBestLikelihood(), "SECOND_LIKELIHOOD="+result.getSecondBestLikelihood(),
				"LIKELIHOOD_DELTA=" + result.getLikelihoodDelta(), "NORMALIZED_LIKELIHOOD="+result.getNormalizedLikelihood(), "SHANNON_WEAVER_DIVERSITY="+divFormat.format(diversity), "SHANNON_WEAVER_EQUITABILITY="+divFormat.format(equitability),
				"KNOWN_DONOR_TAG="+ kdt};

		String h = StringUtils.join(header, "\t");

		try {
			outTerse.write("#");  // comment line for R code.
			outTerse.write(h); outTerse.newLine();
			outTerse.flush();
		} catch (IOException e) {
			throw new RuntimeException("Trouble writing to output file");
		}

	}

	private void writeDonorRatioOutputHeader (final BufferedWriter outTerse) {
		List<String> paths= this.INPUT_BAM.stream().map(x -> x.getAbsolutePath()).collect(Collectors.toList());		
		String bamList=StringUtils.join(paths, ",");
		
		String [] header={"INPUT_BAM="+bamList, "INPUT_VCF="+this.INPUT_VCF.getAbsolutePath(),
				"GQ_THRESHOLD="+Integer.toString(this.GQ_THRESHOLD), "FRACTION_SAMPLES_PASSING="+Double.toString(this.FRACTION_SAMPLES_PASSING),
				"MIN_NUM_VARIANT_SAMPLES="+Integer.toString(this.MIN_NUM_VARIANT_SAMPLES),
				"READ_MQ="+Integer.toString(this.READ_MQ)};
		String h = StringUtils.join(header, "\t");

		try {
			outTerse.write("#");  // comment line for R code.
			outTerse.write(h); outTerse.newLine();
			outTerse.flush();
		} catch (IOException e) {
			throw new RuntimeException("Trouble writing to output file");
		}
	}

	private void writeCommonDonorRatioOutput (final Map<String, Double> summary, final ObjectCounter<String> knownDonorCounts, final CommonSNPsData data, final BufferedWriter outTerse) {
		// DecimalFormat pctFormat = new DecimalFormat("0.###");
		DecimalFormat fracFormat = new DecimalFormat("0.#####");
		if (outTerse==null) return;
		try {
			List<String> header = new ArrayList<>();
			header.addAll(Arrays.asList("DONOR", "REPRESENTATION"));
			if (knownDonorCounts!=null)
				header.addAll(Arrays.asList("KNOWN", "COUNT"));
			if (this.REPORT_ALLELE_COUNTS)
				header.addAll(Arrays.asList("HOM_REF", "HET", "HOM_VAR", "MISSING"));

			Map<String, Double> knownFracMap = convertCountsToFractionalRep(knownDonorCounts);

			String h = StringUtils.join(header, "\t");
			outTerse.write(h); outTerse.newLine();
			List<String> samples = new ArrayList<>(summary.keySet());
			Collections.sort(samples);
			for (String sample: samples) {
				double ratio=summary.get(sample);
				String frac = fracFormat.format(ratio);
				List<String> line = new ArrayList<>();
				line.addAll(Arrays.asList(sample, frac));
				if (knownDonorCounts!=null) {
					Double knownFrac = knownFracMap.get(sample);
					if (knownFrac==null) {
						line.add("0");
						line.add("0");
					}
					else {
						line.add(fracFormat.format(knownFrac));
						line.add(Integer.toString(knownDonorCounts.getCountForKey(sample)));
					}

				}
				if (this.REPORT_ALLELE_COUNTS) {
					ObjectCounter<GenotypeType> counts = data.getCountGenotypes(sample);
					line.add(Integer.toString(counts.getCountForKey(GenotypeType.HOM_REF)));
					line.add(Integer.toString(counts.getCountForKey(GenotypeType.HET)));
					line.add(Integer.toString(counts.getCountForKey(GenotypeType.HOM_VAR)));
					line.add(Integer.toString(counts.getCountForKey(GenotypeType.NO_CALL)));
				}
				String l = StringUtils.join(line, "\t");
				outTerse.write(l); outTerse.newLine();
			}
			outTerse.close();
		} catch (IOException e) {
			throw new RuntimeException("Trouble writing to summary file");
		}
		return;
	}

	private Map<String, Double> convertCountsToFractionalRep (final ObjectCounter<String> counter) {
		if (counter==null) return null;
		Map<String, Double> result = new HashMap<>();
		int total = counter.getTotalCount();
		for (String k: counter.getKeys()) {
			double frac = (double) counter.getCountForKey(k)/(double) total;
			result.put(k, frac);
		}
		return result;
	}
		
	/**
	 * Build the data set of common SNPs for optimization.
	 * @param data The current data set to add a variant to.
	 * @param vc The variant context from the VCF
	 * @param pileup The read data from the BAM.
	 * @return The updated
	 */
	private CommonSNPsData processPileUp (final CommonSNPsData data, final VariantContext vc, final SNPGenomicBasePileUp pileup, final VariantContextWriter vcfWriter) {
		char refBase = StringUtil.byteToChar(vc.getReference().getBases()[0]);
		char altBase = CensusSeqUtils.getAltBase(vc);

		int readsRef = pileup.getCountBase(refBase);
		int readsAlt = pileup.getCountBase(altBase);
		GenotypesContext gc = vc.getGenotypes();
		Iterator<Genotype> iter = gc.iterator();
		int numSamples = gc.size();
		Interval i = pileup.getSNPInterval();

		String [] sampleNames = new String [numSamples];
		int [] altAlleleCounts = new int [numSamples];
		int counter=0;

		if (iter.hasNext()==false)
			log.info("STOP");

		while (iter.hasNext()) {
			Genotype g = iter.next();
			int altAlleleCount=0;
			if (g.isNoCall()) altAlleleCount=-1;
			if (g.isHomVar()) altAlleleCount=2;
			if (g.isHet()) altAlleleCount=1;
			String donorName = g.getSampleName();
			sampleNames[counter]=donorName;
			altAlleleCounts[counter]=altAlleleCount;
			counter++;
		}
		if (counter < numSamples-1)
			log.info("STOP for missing data");
		data.addSNP(sampleNames, altAlleleCounts, readsRef, readsAlt);
		if (vcfWriter!=null) vcfWriter.add(vc);
		return (data);
	}
		
	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new CensusSeq().instanceMain(args));
	}
}
