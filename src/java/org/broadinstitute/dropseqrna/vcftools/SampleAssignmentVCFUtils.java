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
package org.broadinstitute.dropseqrna.vcftools;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.GatherDigitalAlleleCounts;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPInfoCollection;
import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.sampleassignment.AssignCellsToSamples;
import org.broadinstitute.dropseqrna.utils.TransformingIterator;
import org.broadinstitute.dropseqrna.utils.VariantContextProgressLoggerIterator;
import org.broadinstitute.dropseqrna.vcftools.filters.CallRateVariantContextFilter;
import org.broadinstitute.dropseqrna.vcftools.filters.ChromosomeVariantFilter;
import org.broadinstitute.dropseqrna.vcftools.filters.FlipSNPFilter;
import org.broadinstitute.dropseqrna.vcftools.filters.SimpleDiploidVariantContextFilter;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;

/**
 * One stop shopping for consistently getting data from VCFs.
 * @author nemesh
 *
 */
public class SampleAssignmentVCFUtils {
	
	
	/**
	 * Since we iterate through the VCF once to establish an interval list of sites to query, and then once to look at genotypes,
	 * have a 1 stop shop to set up the filtered iterator that handles all processing of VariantContext data.
	 * We require at least one sample to have a called genotype for the variant to enter the iterator.
	 * @param VCFFile The VCF file to read.
	 * @param samples A set of samples to filter on.  Can be empty to not filter samples.
	 * @param log Should a progress logger iterator be added to this iterator chain at the start?  Pass in a Log object if you want logging.
	 * @return A VCF iterator that can (optionally) reduce the number of samples in the VCF, find high quality variants, etc.
     *
	 */
	public static PeekableIterator<VariantContext> getVCFIterator (final VCFFileReader vcfReader, final List<String> samples, final boolean retainMonomorphicSNPs, final int genotypeGCThreshold, final double fractionPassing,  final List<String> ignoredChromosomes, final Log log) {
		Iterator<VariantContext> filteredIter=vcfReader.iterator(); 
		List<String> finalSamples = samples;
		if (samples.size()>0) {
			finalSamples = validateSampleNamesInVCF(vcfReader, samples, log);
			if (sampleListMatchesVcf(vcfReader, finalSamples)) finalSamples=Collections.EMPTY_LIST;
		}
		
		// filter by samples if requested AND if the samples are not the default VCF samples.  Otherwise make samples an empty list.
		
		final PeekableIterator<VariantContext> vcfIterator = filterVCFIterator(filteredIter, finalSamples, retainMonomorphicSNPs, genotypeGCThreshold, fractionPassing, ignoredChromosomes, log);
		return vcfIterator;
	}

	public static PeekableIterator<VariantContext> filterVCFIterator (final Iterator<VariantContext> iter, final List<String> samples, final boolean retainMonomorphicSNPs, final int genotypeGCThreshold, final double fractionPassing,  final List<String> ignoredChromosomes, final Log log) {
		Iterator<VariantContext> filteredIter=iter;
		if (log!=null)
			filteredIter = new VariantContextProgressLoggerIterator(filteredIter, new ProgressLogger(log));
		
		// filter on chromosomes.
		filteredIter = new PeekableIterator<>(new ChromosomeVariantFilter(filteredIter, ignoredChromosomes));

		if (genotypeGCThreshold == -1) {
			log.info("Genotype Quality Filter disabled.  Enabling A/T, C/G SNP Filter to eliminate potential allele flipping variants");
			filteredIter = new FlipSNPFilter(filteredIter);
		}

		if (samples.size()>0)
			filteredIter = new VariantContextSampleTransformer (filteredIter, new HashSet<>(samples));
		
		// filter on the variant #alleles, quality, etc.
		filteredIter = new SimpleDiploidVariantContextFilter(filteredIter, true, true, 2, retainMonomorphicSNPs);
		// filter on genotype call rate
		filteredIter = new CallRateVariantContextFilter(filteredIter, genotypeGCThreshold, fractionPassing);
		
		final PeekableIterator<VariantContext> vcfIterator = new PeekableIterator<>(filteredIter);
		return vcfIterator;
	}
	
	/**
	 * Are the samples in list provided the complete list of samples in the vcfReader, and in the same order?
	 * @param vcfReader
	 * @param samples
	 * @return
	 */
	public static boolean sampleListMatchesVcf (final VCFFileReader vcfReader, final List<String> samples) {
		List<String> vcfSamples=vcfReader.getFileHeader().getSampleNamesInOrder();
		return (vcfSamples.equals(samples));
		
	}

	/**
	 * Scan the VCF, and build a set of SNP intervals to look at in the BAM.
	 * This filters out variants not flagged as "PASSED", variants that aren't SNPs, or variants with more than 2 alleles.
	 * The intervalList generated uses the sequence dictionary from the VCF file.
	 * This is useful as it enforces how the BAM chromosomes are later sorted when the SampleGenotypeProbabilitiesIterator is created.
	 * @param VCFFile The VCF File to parse
	 * @param bamFile The BAM file to get a sequence dictionary for the IntervalList from.
	 * @param preserveIntervalNames Set to true to use the genotype site ID (via site.getID) as the interval name
	 * @param vcfWriter (Optional) write the VCF records from the iterator to this file.  Set to null to ignore.
	 * @return An intervalList of variant locations to examine in the BAM file.
	 */
	public static SNPInfoCollection getSNPInfoCollection (final File vcfFile, final List<String> samples, final boolean retainMonomorphicSNPs, final int genotypeGCThreshold, 
			final double fractionPassing, final List<String> ignoredChromosomes, final VariantContextWriter vcfWriter, final boolean preserveIntervalNames, final Log log) {
		
		final VCFFileReader vcfReader = new VCFFileReader(vcfFile, false);
		// List<String> finalSamples=validateSampleNamesInVCF(vcfReader, samples, log);
		
		PeekableIterator<VariantContext> vcfIterator = getVCFIterator(vcfReader, samples, retainMonomorphicSNPs, genotypeGCThreshold, fractionPassing, ignoredChromosomes, log);
				
		// Begin optional work to remove duplicates		
		SNPInfoCollection result = new SNPInfoCollection(vcfIterator, vcfReader.getFileHeader().getSequenceDictionary(), preserveIntervalNames, vcfWriter, log);				
		return result;
	}

	/**
	 * For an iterator of VCF Records, scan the records and build a set of SNP intervals to look at in the BAM.
	 * The intervalList generated uses the sequence dictionary passed in.
	 * @param vcfIterator The iterator of VCF records
	 * @param sd The sequence dictionary
	 * @param log An optional log file to log to.  Can be null.
	 * @return An intervalList of variant locations from the vcfIterator
	 */
	/*
	public static IntervalList getSNPIntervals  (final Iterator<VariantContext> vcfIterator, final SAMSequenceDictionary sd, final Log log) {
		SAMFileHeader h = new SAMFileHeader();
		h.setSequenceDictionary(sd);
		IntervalList result = new IntervalList(h);
		if (log!=null) log.info("Scanning VCF to find potential SNP sites");
		while (vcfIterator.hasNext()) {
			VariantContext site = vcfIterator.next();
			Interval variant = new Interval(site.getContig(), site.getStart(),site.getEnd(), true, site.getID());
			result.add(variant);
		}
		if (log!=null) log.info("Found [" + result.getIntervals().size() +"] potential SNP sites to query.");
		return (result);
	}
	*/
	/**
	 * Also writes the sites that pass variouts filters in the vcfIterator to be written to a VCF Writer.
	 * @param vcfIterator
	 * @param sd
	 * @param log
	 * @param vcfWriter
	 * @return
	 */
	/*
	public static IntervalList getSNPIntervals  (final Iterator<VariantContext> vcfIterator, final SAMSequenceDictionary sd, final Log log, 
												 final VariantContextWriter vcfWriter, final boolean preserveIntervalNames) {
		SAMFileHeader h = new SAMFileHeader();
		h.setSequenceDictionary(sd);
		IntervalList result = new IntervalList(h);
		if (log!=null) log.info("Scanning VCF to find potential SNP sites");
		while (vcfIterator.hasNext()) {
			VariantContext site = vcfIterator.next();

			Interval variant;
			if (preserveIntervalNames)
				variant = new Interval(site.getContig(), site.getStart(),site.getEnd(), true, site.getID());
			else
				variant = new Interval(site.getContig(), site.getStart(),site.getEnd());
			result.add(variant);
			vcfWriter.add(site);
		}
		if (log!=null) log.info("Found [" + result.getIntervals().size() +"] potential SNP sites to query.");
		return (result);
	}
	*/
	
	/*
	public static IntervalAndFrequencyResult getIntervalAndFrequency (final Iterator<VariantContext> vcfIterator, final SAMSequenceDictionary sd, final Log log, final VariantContextWriter vcfWriter, final boolean preserveIntervalNames) {
		IntervalAndFrequencyResult fullResult = new IntervalAndFrequencyResult();

		SAMFileHeader h = new SAMFileHeader();
		h.setSequenceDictionary(sd);
		IntervalList result = new IntervalList(h);
		if (log!=null) log.info("Scanning VCF to find potential SNP sites");
		while (vcfIterator.hasNext()) {
			VariantContext site = vcfIterator.next();
			Interval variant;
			if (preserveIntervalNames)
				variant = new Interval(site.getContig(), site.getStart(),site.getEnd(), true, site.getID());
			else
				variant = new Interval(site.getContig(), site.getStart(),site.getEnd());
			result.add(variant);
			vcfWriter.add(site);
		}
		if (log!=null) log.info("Found [" + result.getIntervals().size() +"] potential SNP sites to query.");
		return (fullResult);

	}
	*/


	/**
	 * If the set of samples is non-empty, validate how many samples overlap in the VCF file and the set of samples.
	 * If the intersect set is of size 0, something is very wrong so throw an error.
	 * If the intersect doesn't find all of the requested samples, something is very wrong so throw an error.
	 * @param vcfReader The VCF reader
	 * @param samples A set of strings of sample names expected to be in the VCF.  If this is empty, this operation is a no-op.
	 */
	public static List<String> validateSampleNamesInVCF (final VCFFileReader vcfReader, final List<String> samples, final Log log) {
		// short circuit for no-op.
		if (samples.isEmpty()) return samples;
		Set<String> samplesInVCF = new HashSet<>(vcfReader.getFileHeader().getSampleNameToOffset().keySet());
		samplesInVCF.retainAll(samples);
		if (samplesInVCF.size()==0) {
			samplesInVCF = vcfReader.getFileHeader().getSampleNameToOffset().keySet();
			throw new IllegalArgumentException("Didn't find any of these samples from the VCF:"+samplesInVCF +" in the submitted a list of samples " + samples.toString());
		}
		
		if (log!=null) {
			String msg = "Found " + Integer.toString(samplesInVCF.size()) + " samples in VCF and requested sample list out of " + Integer.toString(samples.size()) + " requested";
			log.info(msg);
		}
				
		if (samplesInVCF.size()<samples.size())
			throw new IllegalArgumentException("Did not find all of the requested samples.  Can not continue.");

		List<String> finalSamples = new ArrayList<>();
		for (String s: samples)
			if (samplesInVCF.contains(s)) finalSamples.add(s);
		return (finalSamples);
	}

	/**
	 * Runs subContextFromSamples on each VariantContext record to filter a VariantContext record
	 * to only contain data for the samples supplied to this iterator.
	 * @author nemesh
	 *
	 */
	public static class VariantContextSampleTransformer extends TransformingIterator<VariantContext, VariantContext> {
		private final Set<String> sampleNames;
		public VariantContextSampleTransformer(final Iterator<VariantContext> underlyingIterator, final Set<String> sampleNames) {
			super(underlyingIterator);
			this.sampleNames=sampleNames;
		}

		@Override
		public VariantContext next() {
			VariantContext vc = this.underlyingIterator.next();
			return vc.subContextFromSamples(sampleNames);
		}
	}

	/**
	 * If there's a file with a list of samples to use, parse that and use that subset of samples
	 * Otherwise, get the list of samples in the VCF header.
	 * @param vcfReader
	 * @return
	 */
	public static List <String> getVCFSamples (final VCFFileReader vcfReader, final File vcfSampleFile) {
		List<String> samples;
		if (vcfSampleFile!=null) {
			IOUtil.assertFileIsReadable(vcfSampleFile);
			samples=ParseBarcodeFile.readCellBarcodeFile(vcfSampleFile);
		} else
			samples=vcfReader.getFileHeader().getSampleNamesInOrder();
		return samples;
	}

	/**
	 * Returns a writer, or null if the out file is null.
	 * @param vcfReader
	 * @param out
	 * @return
	 */
	public static VariantContextWriter getVCFWriter(final VCFFileReader vcfReader, final File out) {
		if (out == null)
			return null;
		VariantContextWriterBuilder builder = new VariantContextWriterBuilder()
				.setReferenceDictionary(
						vcfReader.getFileHeader().getSequenceDictionary())
				.setOption(Options.INDEX_ON_THE_FLY).setBuffer(8192)
				.setOutputFile(out);
		// don't need this
		// TODO: remove after testing.
		/*
		if (out.getName().endsWith(".gz"))
			builder.setOutputFileType(OutputType.BLOCK_COMPRESSED_VCF);
		else
			builder.setOutputFileType(OutputType.VCF);
		*/
		VariantContextWriter writer = builder.build();

		return writer;
	}



}
