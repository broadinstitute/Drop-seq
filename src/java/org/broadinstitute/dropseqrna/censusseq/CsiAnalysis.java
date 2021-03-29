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
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.apache.commons.io.FileUtils;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.LikelihoodUtils;
import org.broadinstitute.dropseqrna.censusseq.VCFPileupJointIterator.JointResult;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.VCFUtils;
import org.broadinstitute.dropseqrna.utils.readiterators.SamFileMergeUtil;
import org.broadinstitute.dropseqrna.utils.readiterators.SamHeaderAndIterator;
import org.broadinstitute.dropseqrna.vcftools.SampleAssignmentVCFUtils;
import org.broadinstitute.dropseqrna.vcftools.filters.FindMonomorphicSitesInDonorPool;
import org.broadinstitute.dropseqrna.utils.AssertSequenceDictionaryIntersection;
import org.broadinstitute.dropseqrna.utils.FileListParsingUtils;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.metrics.StringHeader;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

/**
 * Find Fortuitously Absent Common SNPs.
 *
 * @author nemesh 
 *
 */

@CommandLineProgramProperties(summary = "Given a VCF file, a list of donors from that file, and BAM file, emit the number of ref/alt reads for variants that are entirely ref in the pool, and calculate the estimated contamination rate of the unexpected donors in the pool.", oneLineSummary = "Find contamination signals in census-seq data", programGroup = DropSeq.class)
public class CsiAnalysis extends CommandLineProgram {

	private static final Log log = Log.getInstance(CsiAnalysis.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze. This argument can accept wildcards, or a file with the suffix .bam_list that contains the locations of multiple BAM files", minElements = 1)
	public List<File> INPUT_BAM;
	
	@Argument(doc = "The input VCF.")
	public File INPUT_VCF;

	@Argument(doc = "A file containing a list of samples to filter the VCF by.  Has 1 column, 1 entry per row.  Each entry is a single sample name from the VCF.", optional = false)
	public File SAMPLE_FILE;

	@Argument (doc = "An INFO tag variants that indicates the allele frequency.  If this is specified, use this expected allele frequency for the variant "
			+ "instead of generating one based on the donors in the VCF not in the SAMPLE_FILE.", optional=true)
	public String ALLELE_FREQ_TAG;
	
	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output file containing the summarized contamination rate of the experiment.", optional=false)
	public File OUTPUT;

	@Argument(doc = "Output of each SNP record in the data set that passes filters at the VCF level and has at least 1 read in the BAM", optional = true)
	private File OUTPUT_VERBOSE;

	@Argument(doc = "Output the list of sites that were monomprhic in all individuals in the sample file, and pass other filter criteria.", optional=true)
	public File VCF_OUTPUT;

	@Argument(doc = "The minimum minor allele frequency in the general population for a SNP to be considered in the summary.")
	public Double MINIMUM_MAF = 0.025;

	@Argument(doc = "The minimum genotype quality for a variant.  If VCF file has no genotype quality, this will be set to -1 to disable filtering.")
	public Integer GQ_THRESHOLD = 30;

	@Argument(doc = "At least <FRACTION_SAMPLES_PASSING> samples must have genotype scores >= GQ_THRESHOLD for the variant in the VCF to be included in the analysis.")
	public double FRACTION_SAMPLES_PASSING = 0.9;

	@SuppressWarnings({ "unchecked", "rawtypes" })
	@Argument(doc = "A list of chromosomes to omit from the analysis.  The default is to omit the sex chromosomes.")
	public List<String> IGNORED_CHROMOSOMES = new ArrayList(Arrays.asList("X", "Y", "MT"));

	@Argument(doc = "The map quality of the read to be included.")
	public Integer READ_MQ = 10;

	@Argument(doc = "The minimum base quality of the SNP base on the read", optional = true)
	public Integer MIN_BASE_QUALITY = 10;

	@Argument(doc = "The estimated sequencing error rate")
	public double SEQUENCING_ERROR_RATE = 0.001;

	@Argument(doc = "This option is intended to make evaluation of the algorithm easier when supplied proper training data.  "
			+ "When tagging reads from known donors with the KNOWN_DONOR_TAG and synthetically mixing them to evalutate the algorithm, "
			+ "this counts the number of reads observed with the tag per donor.  Non-informative reads (low map quality, PCR duplicates) are not counted. "
			+ "Reads are only counted on contigs for which there is at least one SNP present in the VCF (minus IGNORED_CHROMOSOES that are ignored) "
			+ "that could inform the data. This tends to exclude bait contigs and alternative haplotypes."
			+ "The final contamination rate is based on the number of reads that originate from either reads that are not labeled with this tag, plus "
			+ "any reads that come from donors in the EXCLUDE_KNOWN_DONOR list."
			+ "For example, if there were 40 donors in the BAM file and all reads were tagged with their donor of origin via this tag, and the "
			+ "EXCLUDE_KNOWN_DONOR was set to DONOR_A, the fraction of reads from DONOR_A in the input should be estimated by the output contaimination rate.", optional=true)
	public String KNOWN_DONOR_TAG=null;
	
	@Argument (doc="Exclude one or more KNOWN donors from the calculation of contamination rate.  This removes those donors from the SAMPLE_FILE, and counts the fraction of reads these"
			+ "donors contain in the data set as the contamination rate.", optional=true)
	public List<String> EXCLUDE_KNOWN_DONOR;
	
	private final String SNP_TAG = "YS";
	private final double MISSING_AF_VALUE=-1d;

	DecimalFormat pctFormat = new DecimalFormat("0.####");

	@Override
	public int doWork() {

		if (this.EXCLUDE_KNOWN_DONOR!=null && this.EXCLUDE_KNOWN_DONOR.size()>0 && this.KNOWN_DONOR_TAG==null) {
			log.error("If excluding 1 or more donors, there must be a tag on all reads to indicate the donor of origin of the data!");
			return(1);
		}
		this.INPUT_BAM = FileListParsingUtils.expandFileList(INPUT_BAM);		
		IOUtil.assertFileIsReadable(INPUT_VCF);
		IOUtil.assertFileIsReadable(SAMPLE_FILE);
		IOUtil.assertFileIsWritable(this.OUTPUT);
		if (this.VCF_OUTPUT!=null) IOUtil.assertFileIsWritable(this.VCF_OUTPUT);

		
		// set up the optional output
		BufferedWriter outVerbose = null;
		if (OUTPUT_VERBOSE != null) {
			IOUtil.assertFileIsWritable(OUTPUT_VERBOSE);
			outVerbose = new BufferedWriter(new OutputStreamWriter(IOUtil.openFileForWriting(OUTPUT_VERBOSE)));
		}

		VCFFileReader vcfReader = new VCFFileReader(this.INPUT_VCF, true);
		
		SamHeaderAndIterator headerAndIter= SamFileMergeUtil.mergeInputs(this.INPUT_BAM, false, SamReaderFactory.makeDefault());				
		AssertSequenceDictionaryIntersection.assertIntersectionObjectVcf(headerAndIter.header, "BAM INPUT(S)", this.INPUT_VCF, log);		
		SAMSequenceDictionary sd = vcfReader.getFileHeader().getSequenceDictionary();

		// validate VCF is indexed.
		if (!VCFUtils.hasIndex(this.INPUT_VCF)) return 1;

		// need to establish this state for multiple VCF parsing operations.
		if (!VCFUtils.GQInHeader(vcfReader)) {
			this.GQ_THRESHOLD = -1;
			log.info("Genotype Quality [GQ] not found in header.  Disabling GQ_THRESHOLD parameter");
		}

		List<String> donorsInPool = ParseBarcodeFile.readCellBarcodeFile(SAMPLE_FILE);
		// validate list.  If not valid, exit.
		boolean validList = validateSampleList (vcfReader, donorsInPool);
		if (!validList)
			return 1;

		// remove requested known donors from calculation.
		if (this.EXCLUDE_KNOWN_DONOR!=null && EXCLUDE_KNOWN_DONOR.size()>0) {
			log.info("Excluding known donors from analysis " + this.EXCLUDE_KNOWN_DONOR.toString());
			donorsInPool.removeAll(this.EXCLUDE_KNOWN_DONOR);
			log.info("Remaining Donors in Pool [" +donorsInPool.size() +"]");
		}		
		
		// all of this work has been abstracted into a utility class so it can be reusable.
		FindMonomorphicSitesInDonorPool fmsdp = new FindMonomorphicSitesInDonorPool(donorsInPool, vcfReader, this.GQ_THRESHOLD, this.TMP_DIR.get(0),
				this.ALLELE_FREQ_TAG, this.MINIMUM_MAF, this.IGNORED_CHROMOSOMES, this.FRACTION_SAMPLES_PASSING);

		IntervalList snpIntervals = fmsdp.getIntervalList();
		File outTempVCF = fmsdp.getMinimalVCFFile();
		copyTempVCFToOutput(this.VCF_OUTPUT, outTempVCF);

		vcfReader.close();

		// build the BAM iterator.
		log.info("Finding SNPs in BAM.");		
		SNPGenomicBasePileupIterator pileUpIter = new SNPGenomicBasePileupIterator(headerAndIter, snpIntervals, SNP_TAG,
				READ_MQ, this.IGNORED_CHROMOSOMES, KNOWN_DONOR_TAG, this.MIN_BASE_QUALITY);

		// reset the iterator for use in the full data set. Use the cleaned up
		// set of variants, which should be smaller and faster.
		// technically this probably doesn't need to get pushed through filters again...
		vcfReader = new VCFFileReader(outTempVCF, false);
		PeekableIterator<VariantContext> vcfIterator = fmsdp.getVCFIterator(vcfReader.iterator(), donorsInPool, this.ALLELE_FREQ_TAG,
		 		this.MINIMUM_MAF, this.GQ_THRESHOLD, this.IGNORED_CHROMOSOMES, this.FRACTION_SAMPLES_PASSING);

		
		// walk through the VCF/BAM looking for SNPs and scoring them.
		Set<String> donorsInPoolSet = new HashSet<>(donorsInPool);

		findErrorAlleles(this.MINIMUM_MAF, this.ALLELE_FREQ_TAG, donorsInPoolSet, pileUpIter, vcfIterator, sd, outVerbose,
				this.OUTPUT);

		vcfReader.close();
		return 0;
	}

	private boolean validateSampleList (final VCFFileReader vcfReader, final List<String> donorsInPool) {
		List<String> copy = new ArrayList<>(donorsInPool);
		List<String> validVcfSamples = SampleAssignmentVCFUtils.validateSampleNamesInVCF(vcfReader, donorsInPool, log);
		copy.removeAll(validVcfSamples);
		if (copy.size()>0) {
			log.info("Samples found in sample list but not VCF " + copy.toString()+"");
			return false;
		}
		return true;
	}

	/**
	 * Iterate over the SNP pileups and VCF file. Where you encounter a SNP in
	 * both, process it.
	 *
	 * @param vcfIterator
	 * @param sd
	 */
	private void findErrorAlleles(final double mafThreshold, final String alleleFreqTag, final Set<String> donorsInPool,
			final SNPGenomicBasePileupIterator pileUpIter,
			final PeekableIterator<VariantContext> vcfIterator, final SAMSequenceDictionary sd,
			final BufferedWriter outVerbose, final File out) {

		PeekableIterator<SNPGenomicBasePileUp> peekablePileUpIter = new PeekableIterator<>(pileUpIter);

		CsiMetrics summary = new CsiMetrics();
		summary.SEQUENCING_ERROR_RATE = this.SEQUENCING_ERROR_RATE;

		if (outVerbose!=null) writeVerboseHeader(outVerbose);
		VCFPileupJointIterator jointIter = new VCFPileupJointIterator(peekablePileUpIter, vcfIterator, sd);
		while (jointIter.hasNext()) {
			JointResult jr = jointIter.next();
			SNPGenomicBasePileUp pileup = jr.getPileup();
			summary = processPileUp(mafThreshold, alleleFreqTag, summary, donorsInPool, jr.getVc(), pileup, outVerbose);
		}

		if (outVerbose!=null) CloserUtil.close(outVerbose);

		jointIter.close();
		JointIteratorCounter counter = jointIter.getCounter();

		if (counter.SAMPLE_GENE_ITER != counter.BOTH)
			log.info("Expected counts don't match: Pileup [" + counter.SAMPLE_GENE_ITER + "] both counter ["
					+ counter.BOTH + "]");

		log.info("Evaluted SNPs [" + counter.BOTH + "]");
		log.info("DONE");

		ObjectCounter<String> knownDonorCounts = pileUpIter.getKnownDonorCounts();
		Double knownRate = getKnownRate(knownDonorCounts, donorsInPool);		

		DecimalFormat divFormat = new DecimalFormat("0.####");

		// write summary
		MetricsFile<CsiMetrics, String> file = new MetricsFile<>();
		// add header.
		List<String> h = new ArrayList<>();
		h.add("INPUT_BAM="+ this.INPUT_BAM);
		h.add("INPUT_VCF="+ this.INPUT_VCF);
		h.add("SAMPLE_FILE="+ this.SAMPLE_FILE);
		h.add("MINIMUM_MAF="+ this.MINIMUM_MAF);
		h.add("READ_MQ="+ this.READ_MQ);
		h.add("MIN_BASE_QUALITY="+ this.MIN_BASE_QUALITY);
		if (this.ALLELE_FREQ_TAG!=null)
			h.add("ALLELE_FREQ_TAG="+ this.ALLELE_FREQ_TAG);

		file.addHeader(new StringHeader(StringUtils.join(h, "\t")));

		summary.calculateStats();
		if (knownDonorCounts==null || knownDonorCounts.getSize()==donorsInPool.size()) {
			summary.KNOWN_CONTAMINATION_RATE="NA";
		} else {
			summary.KNOWN_CONTAMINATION_RATE=divFormat.format(knownRate);
		}
		
		file.addMetric(summary);
		file.write(out);
		log.info("Result:\n" + summary.toString());
	}

	/**
	 * Calculate the known contamination rate based on the number of reads detectable in the donors in the pool compared to the total number of donors.
	 * rate = (1- (totalInPool/total))*100;
	 * @param knownDonorCounts
	 * @param donorsInPool
	 * @return The known contamination rate, or null if knownDonorCounts is null.
	 */
	private Double getKnownRate (final ObjectCounter<String> knownDonorCounts, final Set<String> donorsInPool) {
		if (knownDonorCounts==null) return null;

		int totalInPool = 0;
		int total = knownDonorCounts.getTotalCount();
		for (String d: donorsInPool) {
			int count = knownDonorCounts.getCountForKey(d);
			totalInPool+=count;
		}
		double rate = (1- ((double)totalInPool/(double)total))*100;
		return (rate);
	}

	private void writeVerboseHeader(final BufferedWriter outVerbose) {
		String[] header = { "chr", "pos", "ref_base", "alt_base", "MAF_pool", "MAF_pop", "num_ref_bases", "num_alt_bases" };
		String h = StringUtils.join(header, "\t");
		try {
			outVerbose.write(h);
			outVerbose.newLine();
			outVerbose.flush();
		} catch (IOException e) {
			throw new RuntimeException("Trouble writing to FACS per-snp output");
		}
	}

	private void copyTempVCFToOutput (final File vcfOutput, final File vcfTempFile) {
		if (vcfOutput==null) return;

		try {
			FileUtils.copyFile(vcfTempFile, vcfOutput);
		} catch (IOException e) {
			log.error("Failed to copy VCF file to output");
			e.printStackTrace();
		}
	}


	/**
	 * Extract the info we need for the report The contig, the position, the ref
	 * base [common in the pool] the alt base, the SNP frequency in the pool,
	 * the SNP frequency in the population, the number of ref bases observed in
	 * data, the number of alt bases in data.
	 *
	 * @param donorsInPool
	 * @param vc
	 * @param pileup
	 * @param outVerbose
	 */
	private CsiMetrics processPileUp(final double mafThreshold, final String alleleFreqTag, final CsiMetrics summary,
			final Set<String> donorsInPool, final VariantContext vc, final SNPGenomicBasePileUp pileup,
			final BufferedWriter outVerbose) {

		String contig = vc.getContig();
		int start = vc.getStart();

		// System.out.println(contig+":"+start);

		Allele refAllele = vc.getReference();
		List<Allele> otherAlleles = vc.getAlternateAlleles();
		// find the best alternate allele.
		Allele altAllele = null;

		if (otherAlleles.size()==1)
			altAllele=otherAlleles.get(0);
		else { // if there's more than 1 "other" allele, try to select it.
			// prefer the most "used" allele. in the remaining donors in the VCF.
			int bestCount = 0;
			for (Allele a : otherAlleles)
				if (a.isCalled()) {
					int c = vc.getCalledChrCount(a);
					if (c > bestCount)
						altAllele = a;
				}
		}
		// The donorPool should be monomorphic. It can be the "ref" allele in
		// the population or the alt allele.
		// we want to orient the ref and alt alleles so the common allele in the
		// pool is the "ref"
		// if the countRefAlleleFirstDonor is greater than 0, then the ref and
		// alt bases are in the right orientation.
		// look through all the donors for a SNP that is called, and get its
		// genotype.
		int countRefAlleleFirstDonor = 0;
		VariantContext donorPoolVC = vc.subContextFromSamples(donorsInPool);
		for (String d : donorsInPool) {
			Genotype g = donorPoolVC.getGenotype(d);
			if (g.isHom() || g.isHet()) {
				countRefAlleleFirstDonor = g.countAllele(refAllele);
				break; // you found one donor, you're done.
			}
		}

		// swap the ref and alt base if the pool is 100% alt allele.
		if (countRefAlleleFirstDonor == 0) {
			// the 3 way swap.
			Allele temp = refAllele;
			refAllele = altAllele;
			altAllele = temp;
		}

		char poolRefBase = StringUtil.byteToChar(refAllele.getBases()[0]);
		char poolAltBase = StringUtil.byteToChar(altAllele.getBases()[0]);

		int readsRef = pileup.getCountBase(poolRefBase);
		int readsAlt = pileup.getCountBase(poolAltBase);

		// if the only base seen is not the expected ref or alt allele, do NOT
		// emit this result.
		if (readsRef == 0 && readsAlt == 0)
			return summary;

		double afPool = getAlleleFrequency(donorPoolVC, refAllele);
		double afPoolMinor = 1 - afPool;

		double afAllMinor=this.MISSING_AF_VALUE;

		// If the ALLELE_FREQ_TAG is set,
		if (alleleFreqTag!=null) {
			afAllMinor=vc.getAttributeAsDouble(alleleFreqTag, this.MISSING_AF_VALUE);
			if (afAllMinor==this.MISSING_AF_VALUE)
				log.error("A variant record with a missing allele frequency tag passed through filters.  This shouldn't happen. "+ vc.toString());
			// if the reference base was flipped, flip the minor allele frequency.
			if (countRefAlleleFirstDonor==0)
				afAllMinor = 1-afAllMinor;
				// log.info("Flipped allele frequency tag of SNP" + vc.toStringWithoutGenotypes());

		}
		else {
			double afAll = getAlleleFrequency(vc, refAllele);
			afAllMinor = 1 - afAll;
		}

		// update summary if appropriate.
		if (afAllMinor >= mafThreshold) {
			summary.NUM_SNPS++;
			summary.REF_COUNT += readsRef;
			summary.ALT_COUNT += readsAlt;
			// summary.addAlleleFreq(afAllMinor);
			
			List<Byte> bases = pileup.getBases();
			List<Byte> quals = pileup.getQualities();
			
			for (int i=0; i<bases.size(); i++) {
				// add the allele frequency for each base.  This weighs SNPs that have many observations more.
				summary.addAlleleFreq(afAllMinor);
				Byte base = bases.get(i);
				Byte quality= quals.get(i);				
				double errorRate = LikelihoodUtils.getInstance().phredScoreToErrorProbability(quality);				
				summary.addBaseErrorProbability(errorRate);
				// add in the base quality aware info!
				// if the base observed is the reference base, then construct the dosages of ref and alt where alt is the error.
				/*
				if (base==refAllele.getBases()[0]) {
					summary.REF_DOSAGE+=(1-errorRate);
					summary.ALT_DOSAGE+=errorRate;
				} // if the base observed is the reference base, then construct the dosages of ref and alt where ref is the error.
				else if (base==altAllele.getBases()[0]) {
					summary.REF_DOSAGE+=errorRate;
					summary.ALT_DOSAGE+=(1-errorRate);
				}
				*/
				// we ignore the case where the base isn't the ref or the alt allele.  
				//TODO: can we do anything with this info?  It doesn't support the ref base, but we can't use the MAF of the SNP because this is technically another site.
			}						
		}

		String[] line = { contig, Integer.toString(start), Character.toString(poolRefBase),
				Character.toString(poolAltBase), pctFormat.format(afPoolMinor), pctFormat.format(afAllMinor),
				Integer.toString(readsRef), Integer.toString(readsAlt) };

		if (outVerbose!=null)
			try {
				outVerbose.write(StringUtils.join(line, "\t"));
				outVerbose.newLine();
				outVerbose.flush();
			} catch (IOException e) {
				throw new RuntimeException("Trouble writing to FACS per-snp output");
			}
		return summary;
	}

	public static double getAlleleFrequency(final VariantContext vc, final Allele refAllele) {
		int total = 0;
		int refCount = 0;
		GenotypesContext gc = vc.getGenotypes();
		Iterator<Genotype> iter = gc.iterator();
		while (iter.hasNext()) {
			Genotype g = iter.next();
			if (!g.isNoCall()) {
				refCount += g.countAllele(refAllele);
				total += g.getAlleles().size();
			}
		}
		return (double) refCount / (double) total;
	}


	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new CsiAnalysis().instanceMain(args));
	}
}
