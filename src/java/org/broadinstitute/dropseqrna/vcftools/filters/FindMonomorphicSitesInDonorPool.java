/*
 * MIT License
 *
 * Copyright 2018 Broad Institute
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

package org.broadinstitute.dropseqrna.vcftools.filters;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import htsjdk.samtools.util.*;

import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPInfoCollection;
import org.broadinstitute.dropseqrna.censusseq.CensusSeqUtils;
import org.broadinstitute.dropseqrna.utils.VariantContextProgressLoggerIterator;
import org.broadinstitute.dropseqrna.vcftools.SampleAssignmentVCFUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * Abstracts out parsing through a VCF file and BAM file to find sites that are high quality AND monomorphic in the donor pool.
 * Returns both an IntervalList of sites, as well as a temp file that contains the list of sites that are valid in the interval list.
 * @author nemesh
 *
 */
public class FindMonomorphicSitesInDonorPool {

	private static final Log log = Log.getInstance(FindMonomorphicSitesInDonorPool.class);
	private final List<String> donorsInPool;
	private final VCFFileReader vcfReader;
	private final File tmpDir;
	private final int genotypeQuality;
	private final String alleleFrequencyTag;
	private final Double minimum_MAF;
	private final List<String> ignoredChromosomes;
	private final double fractionDonorsPassing;

	// the outputs.
	private IntervalList snpIntervals;
	private File outputVCF;

	public FindMonomorphicSitesInDonorPool (final Collection <String> donorsInPool, final VCFFileReader vcfReader, final int genotypeQuality, final File tmpDir, final String alleleFrequencyTag,
			final Double minimum_MAF, final List<String> ignoredChromosomes, final double fractionDonorsPassing) {
		this.donorsInPool=new ArrayList<>(donorsInPool);
		this.vcfReader=vcfReader;
		this.genotypeQuality=genotypeQuality;
		this.tmpDir=tmpDir;
		this.alleleFrequencyTag=alleleFrequencyTag;
		this.minimum_MAF=minimum_MAF;
		this.ignoredChromosomes=ignoredChromosomes;
		this.fractionDonorsPassing=fractionDonorsPassing;
		this.snpIntervals=null;
		this.outputVCF=null;


	}

	public IntervalList getIntervalList () {
		if (snpIntervals==null)
			getSitesRefInPool();
		return snpIntervals;
	}

	public File getMinimalVCFFile () {
		if (this.outputVCF==null)
			getSitesRefInPool();
		return this.outputVCF;
	}

	/**
	 * Finds a list of SNP intervals that are homozygous (same direction) in all donors in the pool
	 * as defined by the sampleFile.
	 */
	private void getSitesRefInPool () { 

		// set up the variant writer to write to either the output file OR a temp file.
		// write to a temp directory on the first pass so we can iterate on it
		// along with the sorted BAM.

		try {
			this.outputVCF = File.createTempFile("tmp_vcf_", ".txt.gz", tmpDir);
			this.outputVCF.deleteOnExit();
		} catch (IOException e) {
			throw new RuntimeIOException("Exception creating temp file in " + tmpDir.getAbsolutePath(), e);
		}
		log.info("Writing temp VCF to " + outputVCF.getAbsolutePath());
		final VariantContextWriter vcfWriter = CensusSeqUtils.getVCFWriter(outputVCF, vcfReader, new ArrayList<>(vcfReader.getFileHeader().getSampleNameToOffset().keySet()));

		log.info("Looking through VCF for SNPs that fit criteria.  Will search for these in BAM.");
		// find sites in the VCF.
		final PeekableIterator<VariantContext>  vcfIterator = getVCFIterator(vcfReader.iterator(), donorsInPool, alleleFrequencyTag, minimum_MAF, genotypeQuality, ignoredChromosomes, fractionDonorsPassing);

		// last argument is a bit funky. If you aren't using the common SNP
		// analysis, you need to preserve the full interval names.
		final SNPInfoCollection snpIntervals = new SNPInfoCollection(vcfIterator, vcfReader.getFileHeader().getSequenceDictionary(), true, vcfWriter, log);				
		vcfIterator.close();
		vcfWriter.close();
		this.snpIntervals=snpIntervals.getIntervalList();
	}

	/**
	 * Iterate over a VCF iterator with a number of filters in place, only returning valid records in the iterator.
	 * @param iter An existing iterator over the VCF file
	 * @param donorsInPool a list of donors that will be tested.
	 * @param alleleFreqTag a tag on the VCF records to indicate the allele frequency in a reference population.  Can be null.
	 * @param minMAF A minimum minor allele frequency in the population (entire VCF or reference population.)
	 * @param genotypeQuality The minimum genotype quality to consider a genotype.
	 * @param ignoredChromosomes Chromosomes to eliminate from analysis
	 * @param fractionDonorsPassing What fraction of the total number of donors in the population AND what fraction of donors in the pool must both be called with a genotype (not no-call genotype, GQ>= threshold if GQ is available)
	 * @return A new VCF iterator with filters applied.
	 */
	public PeekableIterator<VariantContext> getVCFIterator(final Iterator<VariantContext> iter, final List<String> donorsInPool, final String alleleFreqTag, final Double minMAF, final int genotypeQuality, final List<String> ignoredChromosomes, final double fractionDonorsPassing) {

		Iterator<VariantContext> filteredIter=iter;
		filteredIter = new VariantContextProgressLoggerIterator(filteredIter, new ProgressLogger(log));

		// filter on the variant #alleles, quality, etc.
		filteredIter = new SimpleDiploidVariantContextFilter(filteredIter, true, true, 2);
		// filter on genotype call rate across the entire experiment AND across the pool.
		filteredIter = new CallRateVariantContextFilter(filteredIter, genotypeQuality, fractionDonorsPassing, donorsInPool);
		filteredIter = new CallRateVariantContextFilter(filteredIter, genotypeQuality, fractionDonorsPassing);
		// if the genotype quality is disabled, be extra careful and filter out flip-snps (A/T, C/G SNPs that can be easily messed up in VCFs)
		if (genotypeQuality==-1) {
			log.info("Genotype Quality Filter disabled.  Enabling A/T, C/G SNP Filter to eliminate potential allele flipping variants");
			filteredIter = new FlipSNPFilter(filteredIter);
		}

		// filter problematic chromosomes as defined by the user, usually the sex chromsomes and MT
		filteredIter = new ChromosomeVariantFilter(filteredIter, ignoredChromosomes);

		// retain SNPs monomorphic in the pool of donors
		filteredIter = new MonomorphicOnlyVariantContextFilter(filteredIter, donorsInPool);

		// vcfIterator = new PeekableIterator<>(new ChromosomeVariantFilter(vcfIterator, this.IGNORED_CHROMOSOMES));
		// if using the allele frequency tag, filter out records that don't have the tag set.
		if (alleleFreqTag!=null)
			filteredIter=new AlleleFrequencyTagFilter(filteredIter, alleleFreqTag, minMAF);

		return new PeekableIterator<>(filteredIter);
	}


	/**
	 * This is the original version of the code that has extra hooks for features that were experimental that I may not use.
	 */

	/**
	private void getSitesRefInPoolExperimental () {

		if (!CensusSeqUtils.GQInHeader(vcfReader)) {
			genotypeQuality = -1;
			log.info("Genotype Quality [GQ] not found in header.  Disabling GQ_THRESHOLD parameter");
		}

		SAMSequenceDictionary sd = vcfReader.getFileHeader().getSequenceDictionary();

		// get the list of samples.
		// We want to look at all SNPs that are high quality, regardless of samples.
		// gather up all donors in the VCF.
		// List<String> samplesInVCF = new ArrayList<>(vcfReader.getFileHeader().getSampleNameToOffset().keySet());
		// List<String> samplesInVCF= SampleAssignmentVCFUtils.validateSampleNamesInVCF(vcfReader, donorsInPool, log);

		// set up the variant writer to write to either the output file OR a temp file.
		// write to a temp directory on the first pass so we can iterate on it
		// along with the sorted BAM.
		VariantContextWriter vcfWriter = null;

		try {
			outputVCF = File.createTempFile("tmp_vcf_", ".txt.gz", tmpDir);
			outputVCF.deleteOnExit();
			log.info("Writing temp VCF to " + outputVCF.getAbsolutePath());
			vcfWriter = CensusSeqUtils.getVCFWriter(outputVCF, vcfReader, new ArrayList<>(vcfReader.getFileHeader().getSampleNameToOffset().keySet()));

			if (snpSitesFile==null)
				// you need all the donors in the original VCF to calculate allele frequencies.
				vcfWriter = CensusSeqUtils.getVCFWriter(outTempVCF, vcfReader, new ArrayList<>(vcfReader.getFileHeader().getSampleNameToOffset().keySet()));
			else {
				// open the vcf Site reader just long enough to get header info out.
				VCFFileReader vcfSiteReader = new VCFFileReader(this.SNP_SITES, false);
				vcfWriter = CensusSeqUtils.getVCFWriter(outTempVCF, vcfReader, vcfSiteReader, samplesInVCF);
				vcfSiteReader.close();
			}

		} catch (IOException e) {
			e.printStackTrace();
		}

		PeekableIterator<VariantContext> vcfIterator = null;

		// enhance VCF Iterator with additional records that are assumed to be entirely reference for the donors, and have the allele frequency of the SNP sites.
		// this will produce a Iterator that integrates results from the VCF and SNP site files on the fly.

		SNPSiteEnhancedVCFIterator enhancedIter= null;
		if (this.SNP_SITES!=null) {
			// enhance the VCF file with the sites file.
			enhancedIter = new SNPSiteEnhancedVCFIterator(this.INPUT_VCF, this.SNP_SITES, this.GQ_THRESHOLD, samplesInVCF, this.ANCESTRAL_ALLELE_TAG);
			log.info("Looking through VCF for SNPs that fit criteria.  Will search for these in BAM.");
			vcfIterator = getVCFIterator(enhancedIter, donorsInPool, ALLELE_FREQ_TAG, MINIMUM_MAF);

		} else {
			// open iterator and scan VCF once to get a list of sites.
			log.info("Looking through VCF for SNPs that fit criteria.  Will search for these in BAM.");
			// find sites in the VCF.
			vcfIterator = getVCFIterator(this.INPUT_VCF, donorsInPool, ALLELE_FREQ_TAG, MINIMUM_MAF);
		}

		log.info("Looking through VCF for SNPs that fit criteria.  Will search for these in BAM.");
		// find sites in the VCF.
		vcfIterator = getVCFIterator(vcfReader.iterator(), donorsInPool, alleleFrequencyTag, minimum_MAF, genotypeQuality, ignoredChromosomes, fractionDonorsPassing);

		// last argument is a bit funky. If you aren't using the common SNP
		// analysis, you need to preserve the full interval names.
		final IntervalList snpIntervals = SampleAssignmentVCFUtils.getSNPIntervals(vcfIterator, sd, log, vcfWriter, true);
		vcfIterator.close();
		vcfWriter.close();
		this.snpIntervals=snpIntervals;
	}
	*/

}
