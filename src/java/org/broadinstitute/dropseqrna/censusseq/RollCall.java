package org.broadinstitute.dropseqrna.censusseq;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPInfoCollection;
import org.broadinstitute.dropseqrna.censusseq.VCFPileupJointIterator.JointResult;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.AssertSequenceDictionaryIntersection;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.VCFUtils;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
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
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(summary = "Given a VCF file, a list of donors from that file, and BAM file, find the relative fraction of each donor contained in the sequencing experiment.  Handles both single and paired end reads. By default this filters reads with the PCR Duplicate flag.  ", oneLineSummary = "Determine the donors in a sample pool using private variants", programGroup = DropSeq.class)
public class RollCall extends CommandLineProgram {

	private static final Log log = Log.getInstance(RollCall.class);

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

	@Argument(doc = "The minimum genotype quality for a variant.  Set this value to 0 to not filter by GQ scores if they are present, or to -1 to completely "
			+ "ignore GQ values if they are not set in the genotype info field.  If the GQ field is not set in the VCF header, this will be set to -1 by default.")
	public Integer GQ_THRESHOLD = 30;
	
	@Argument (doc="At least <FRACTION_SAMPLES_PASSING> samples must have genotype scores >= GQ_THRESHOLD for the variant in the VCF to be included in the analysis.")
	public double FRACTION_SAMPLES_PASSING=0.9;

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

	private final String SNP_TAG = "YS";

	@Argument(doc="Output of each SNP record with the sample that is variant and the allele counts.  Only works for private SNP analysis (COMMON_SNP_ANALYSIS=false)", optional=true)
	public File OUTPUT_VERBOSE;
	
	@Argument(doc="If a SNP has an unusual number of reads, skip it.", optional=true)
	public Integer MAX_READS_PER_SNP=null;

	
	private Integer MIN_NUM_VARIANT_SAMPLES=-1;
	
	@Override
	public int doWork() { 
		
		this.INPUT_BAM = FileListParsingUtils.expandFileList(INPUT_BAM);
		IOUtil.assertFileIsReadable(INPUT_VCF);

		// validate VCF is indexed.
		if (!VCFUtils.hasIndex(this.INPUT_VCF)) return 1;

		final VCFFileReader vcfReader = new VCFFileReader(this.INPUT_VCF, true);

		SamHeaderAndIterator headerAndIter= SamFileMergeUtil.mergeInputs(this.INPUT_BAM, false, SamReaderFactory.makeDefault());
		SAMSequenceDictionary sd = vcfReader.getFileHeader().getSequenceDictionary();
		AssertSequenceDictionaryIntersection.assertIntersection(headerAndIter.header, "BAM INPUT(S)", sd, "VCF INPUT", log);		
		
		if (!VCFUtils.GQInHeader(vcfReader)) {
			this.GQ_THRESHOLD=-1;
			log.info("Genotype Quality [GQ] not found in header.  Disabling GQ_THRESHOLD parameter");
		}
		
		// get the list of samples.  Start with the list of all samples, filter to the sample list if needed.  
		List<String> vcfSamples  = CensusSeqUtils.getFinalSamplelist (vcfReader, this.SAMPLE_FILE);
		if (vcfSamples==null) return 1; // exit if the vcf sample list is null
		
		// set up the variant writer to write to either the output file OR a temp file.
		// write to a temp directory on the first pass so we can iterate on it along with the sorted BAM.
		File outTempVCF=CensusSeqUtils.getTempVCFFile (this.TMP_DIR.get(0));
		VariantContextWriter vcfWriter= CensusSeqUtils.getVCFWriter (outTempVCF, vcfReader, vcfSamples);
				
		// set up the optional output
		BufferedWriter outVerbose = null;
		if (OUTPUT_VERBOSE!=null) {
			IOUtil.assertFileIsWritable(OUTPUT_VERBOSE);
			outVerbose = new BufferedWriter(new OutputStreamWriter(IOUtil.openFileForWriting(OUTPUT_VERBOSE)));
			writeVerboseHeader(outVerbose);

		}

		// open iterator and scan VCF once to get a list of sites.
		log.info("Looking through VCF for SNPs that fit criteria.  Will search for these in BAM.");
		// this filters to SNPs in the sample set.
		PeekableIterator<VariantContext> vcfIterator = CensusSeqUtils.getVCFIterator (this.INPUT_VCF, vcfSamples, this.MIN_NUM_VARIANT_SAMPLES, this.GQ_THRESHOLD, 
				this.FRACTION_SAMPLES_PASSING, this.IGNORED_CHROMOSOMES, log, true);
		
		// last argument is a bit funky.  If you aren't using the common SNP analysis, you need to preserve the full interval names.
		// final IntervalList snpIntervals = SampleAssignmentVCFUtils.getSNPIntervals(vcfIterator, sd, log, vcfWriter, true);
		final SNPInfoCollection snpIntervals = new SNPInfoCollection(vcfIterator, vcfReader.getFileHeader().getSequenceDictionary(), true, vcfWriter, log);

		// a requested early exit if there are no SNPs.
		if (snpIntervals.isEmpty()) {
			log.error("No SNP intervals detected!  Check to see if your VCF filter thresholds are too restrictive!");
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
		// this explicitly filters to the best SNP seen on a read if there are multiple SNPs touched by a read
		Map<Interval, Double> genotypeQuality = snpIntervals.getAverageGQ();						
		SNPGenomicBasePileupIterator pileUpIter = new SNPGenomicBasePileupIterator(headerAndIter, snpIntervals.getIntervalList(), SNP_TAG, READ_MQ, this.IGNORED_CHROMOSOMES, KNOWN_DONOR_TAG, genotypeQuality, this.MIN_BASE_QUALITY);

		// reset the iterator for use in the full data set.  Use the cleaned up set of variants, which should be smaller and faster.
		vcfIterator = CensusSeqUtils.getVCFIterator (outTempVCF, vcfSamples, this.MIN_NUM_VARIANT_SAMPLES, this.GQ_THRESHOLD, 
				this.FRACTION_SAMPLES_PASSING, this.IGNORED_CHROMOSOMES, log, true);
				
		// set up the VCF writer for the final pass, if the VCF_OUTPUT is not null.
		vcfWriter = CensusSeqUtils.getVCFWriter (this.VCF_OUTPUT, vcfReader, vcfSamples);

		// run the analysis using only private SNPs.
		int runResult=0;
		
		runResult=findDonorsInBulkRatioPrivateSNPs(vcfSamples, vcfIterator, pileUpIter, snpIntervals.getIntervalList().size(), vcfWriter, sd, outVerbose, outSummary);

		vcfReader.close();
		return runResult;
	}
	
	
	/**
	 * Calculate the donor ratios using private SNPs (only observed on a single individual per SNP.)
	 *
	 * @param vcfIterator The iterator over VCF records
	 * @param vcfWriter A optional writer for informative VCF records 
	 * @param sd The sequence dictionary from the VCF file.
	 * @param outVerbose An optional output writer for verbose results.
	 * @param outSummary 
	 * @return
	 */
	private int findDonorsInBulkRatioPrivateSNPs (final List<String> vcfSamples, final PeekableIterator<VariantContext> vcfIterator, final SNPGenomicBasePileupIterator pileUpIter,
			final int numPossibleSNPs, final VariantContextWriter vcfWriter, final SAMSequenceDictionary sd, final BufferedWriter outVerbose, final BufferedWriter outSummary) {

		PeekableIterator<SNPGenomicBasePileUp> peekablePileUpIter = new PeekableIterator<>(pileUpIter);
		Map<String, SummaryPileUp> summaryResult = new HashMap<>();
		List<SNPSampleRecord> allSNPs = new ArrayList<>();
		Set<String> donorNames = new HashSet<>();

		VCFPileupJointIterator<SNPGenomicBasePileUp> jointIter = new VCFPileupJointIterator<>(peekablePileUpIter, vcfIterator, sd);

		while (jointIter.hasNext()) {
			VCFPileupJointIterator<SNPGenomicBasePileUp>.JointResult jr = jointIter.next();
			SNPSampleRecord r = processPileUp(jr.getVc(), jr.getPileup(), true);
			donorNames.add(r.getDonorName());
			addToSummary(summaryResult, r);
			writeVerboseBody(outVerbose, r);
			if (vcfWriter!=null) vcfWriter.add(jr.getVc());
			allSNPs.add(r);
		}

		jointIter.close();
		JointIteratorCounter counter = jointIter.getCounter();

		if (counter.SAMPLE_GENE_ITER!=counter.BOTH)
			log.info("Expected counts don't match: Pileup [" + counter.SAMPLE_GENE_ITER + "] both counter [" + counter.BOTH+ "]");

		if (vcfWriter!=null) vcfWriter.close();

		if (OUTPUT_VERBOSE!=null) CloserUtil.close(outVerbose);

		Map<String,Double> normalizedRatios = getSampleRatios(summaryResult, true);
		Map<String,Double> rawScores = getSampleRatios(summaryResult, false);
		Map<String, Integer> numSNPs = getNumSNPs(summaryResult);
		ObjectCounter<String> knownResults = pileUpIter.getKnownDonorCounts();

		writePrivateDonorRatioOutputAdditionalData(outSummary, numPossibleSNPs,  counter.BOTH, this.KNOWN_DONOR_TAG);
		writePrivateDonorRatioOutput (vcfSamples, normalizedRatios, rawScores, numSNPs, knownResults, outSummary);

		log.info("Processed [" + counter.BOTH + "] SNPs in BAM + VCF");

		return 0;

	}

	private Map<String, Integer> getNumSNPs (final Map<String, SummaryPileUp> summaryResult) {
		Map<String, Integer> result = new HashMap<>();
		for (String k: summaryResult.keySet()) {
			int v = summaryResult.get(k).getNumSNPs();
			result.put(k, v);
		}
		return result;
	}
	
	/**
	 * Maybe a simpler implementation would be to put the values in the map and keep a running sum of the fractions,
	 * then iterate through again and divide each result by that sum to normalize to 1.
	 * @param summary
	 * @return
	 */
	private Map<String, Double> getSampleRatios (final Map<String, SummaryPileUp> summary, final boolean normalizeToOne) {
		Map<String, Double> result = new HashMap<>();

		List<String> sampleIDs = new ArrayList<>(summary.keySet());
		// hold the ratios before normalizing to 1.
		double [] ratios = new double [summary.size()];
		int counter=0;
		for (String sample: sampleIDs) {
			SummaryPileUp sp=summary.get(sample);
			double ratio = sp.getRatioByAltAlleleFraction();
			ratios[counter]=ratio;
			counter++;
			//result.put(sample, ratio);
		}
		// normalize to 1.
		if (normalizeToOne)
			ratios=OptimizeSampleRatiosGradientFunction.normalizeRatiosToOne(ratios);
		counter=0;
		for (String sample: sampleIDs) {
			result.put(sample, ratios[counter]);
			counter++;
		}
		return result;
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

	void addToSummary (final Map<String, SummaryPileUp> summaryResult, final SNPSampleRecord r) {
		String d = r.getDonorName();
		SummaryPileUp p = summaryResult.get(d);
		if (p==null)
			p = new SummaryPileUp(d);

		p.incrementRefCount(r.getRefCount());
		p.incrementAltCount(r.getAltCount());
		p.incrementNumSNPs();
		summaryResult.put(d, p);
	}

	
	private SNPSampleRecord processPileUp (final VariantContext vc, final SNPGenomicBasePileUp pileup, final boolean hetVarOnly) {
		Interval i = pileup.getSNPInterval();
		List<String> samples = getVariantSamples(vc, hetVarOnly);
		if (samples.size()>1) throw new IllegalStateException("Should only have 1 variant sample per record!");
		String sample = samples.get(0);

		char refBase = StringUtil.byteToChar(vc.getReference().getBases()[0]);
		char altBase = CensusSeqUtils.getAltBase(vc);
		Genotype g = vc.getGenotype(sample);
		
		String strGeno = "NA";
		if (g.isHomRef()) strGeno="ref";
		if (g.isHomVar()) strGeno="alt";
		if (g.isHet()) strGeno = "het";
		SNPSampleRecord result = new SNPSampleRecord(i, refBase, altBase, pileup.getCountBase(refBase),
				pileup.getCountBase(altBase), sample, strGeno, pileup.getQualities());
		return (result);
	}
	

	public List<String> getVariantSamples (final VariantContext vc, final boolean hetVarOnly) {
		GenotypesContext gc = vc.getGenotypes();
		Iterator<Genotype> iter = gc.iterator();
		List<String> result = new ArrayList<>();
		while (iter.hasNext()) {
			Genotype g = iter.next();
			if ( (hetVarOnly & g.isHet()) || !hetVarOnly & !g.isHomRef() )
				result.add(g.getSampleName());
		}
		return (result);
	}

	private void writeVerboseHeader (final BufferedWriter out) {
		if (out==null) return;
		try {
			String [] header = {"CHR", "POS", "REF_BASE", "ALT", "REF_COUNT", "ALT_COUNT", "DONOR", "GENOTYPE", "BASE_QUALITY"};
			String h = StringUtils.join(header, "\t");
			out.write(h);
			out.newLine();
			out.flush();
		} catch (IOException e) {
			throw new RuntimeException("Trouble writing to donor estimation file");
		}
	}

	private void writeVerboseBody (final BufferedWriter out, final SNPSampleRecord r) {
		if (out==null) return;
		try {			
			List<String> header = new ArrayList<String> (Arrays.asList(r.getInterval().getContig(), Integer.toString(r.getInterval().getStart()),
					Character.toString(r.getRefBase()), Character.toString(r.getAltBase()),
					Integer.toString(r.getRefCount()), Integer.toString(r.getAltCount()),
					r.getDonorName(), r.getGenotype()));
			List<Byte> bq = r.getQualities();
			if (bq!=null) 
				header.add(StringUtils.join(bq, ","));
			else
				header.add("");
			String h = StringUtils.join(header, "\t");
			out.write(h);
			out.newLine();
			out.flush();
		} catch (IOException e) {
			throw new RuntimeException("Trouble writing to donor estimation file");
		}
	}

	private void writePrivateDonorRatioOutputAdditionalData (final BufferedWriter outTerse, final int numPossibleSNPs, final int numSNPsFound, final String knownDonorTag) {
		String kdt= knownDonorTag;
		if (kdt==null) kdt= "NULL";

		String [] header={"NUM_POSSIBLE_SNPS="+numPossibleSNPs, "NUM_SNPS_USED="+numSNPsFound, "KNOWN_DONOR_TAG="+ kdt};
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


	private void writePrivateDonorRatioOutput (final List<String> vcfSamples, final Map<String, Double> summary,  final Map<String, Double> rawScores, final Map<String, Integer> numSNPs, final ObjectCounter<String> knownDonorCounts, final BufferedWriter outTerse) {
		// DecimalFormat pctFormat = new DecimalFormat("0.###");
		DecimalFormat fracFormat = new DecimalFormat("0.#####");
		if (outTerse==null) return;
		try {
			List<String> header = new ArrayList<>();
			header.addAll(Arrays.asList("DONOR", "REPRESENTATION", "REP_IRVs", "NUM_SNPS"));
			if (knownDonorCounts!=null)
				header.addAll(Arrays.asList("KNOWN", "COUNT"));

			String h = StringUtils.join(header, "\t");
			outTerse.write(h); outTerse.newLine();
			// List<String> samples = new ArrayList<>(summary.keySet());
			List<String> samples = vcfSamples;
			Collections.sort(samples);
			Map<String, Double> knownFracMap = convertCountsToFractionalRep(knownDonorCounts);

			for (String sample: samples) {
				List<String> line = new ArrayList<>();
				// empty result.
				if (!summary.containsKey(sample))
					line.addAll(Arrays.asList(sample, "0", "0", "0"));
				else {
					Double ratio=summary.get(sample);
					String frac = fracFormat.format(ratio);
					String score = fracFormat.format(rawScores.get(sample));
					int numSNP = numSNPs.get(sample);
					line.addAll(Arrays.asList(sample, frac, score, Integer.toString(numSNP)));
				}
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
				String l = StringUtils.join(line, "\t");
				outTerse.write(l); outTerse.newLine();
			}
			outTerse.close();
		} catch (IOException e) {
			throw new RuntimeException("Trouble writing to summary file");
		}
		return;
	}


	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new CensusSeq().instanceMain(args));
	}
}
