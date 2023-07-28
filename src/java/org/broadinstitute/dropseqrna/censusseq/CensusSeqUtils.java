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

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Set;

import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.utils.IntervalTagComparator;
import org.broadinstitute.dropseqrna.utils.VariantContextSingletonFilter;
import org.broadinstitute.dropseqrna.vcftools.SampleAssignmentVCFUtils;
import org.broadinstitute.dropseqrna.vcftools.filters.ChromosomeVariantFilter;
import org.broadinstitute.dropseqrna.vcftools.filters.CommonVariantContextFilter;
import org.broadinstitute.dropseqrna.vcftools.filters.MonomorphicVariantContextFilter;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

public class CensusSeqUtils {


	private static final Log log = Log.getInstance(CensusSeqUtils.class);
	
	/**
	 * Constructs a VCF writer if the outfile is not null.
	 * New output will only contain the samples listed in vcfSamples.
	 * @param out
	 * @param vcfReader
	 * @param vcfSamples
	 * @return
	 */
	public static VariantContextWriter getVCFWriter (final File out, final VCFFileReader vcfReader, final List<String> vcfSamples) {
		if (out==null) return null;
		IOUtil.assertFileIsWritable(out);
		VariantContextWriter vcfWriter = SampleAssignmentVCFUtils.getVCFWriter(vcfReader, out);

		VCFHeader header = vcfReader.getFileHeader();
		Set<VCFHeaderLine> metaData = header.getMetaDataInInputOrder();
		VCFHeader newHeader = new VCFHeader(metaData, vcfSamples); // set up the new header with the restricted list of samples.
		vcfWriter.writeHeader(newHeader);
		return (vcfWriter);
	}

	public static PeekableIterator<VariantContext> getVCFIterator (final File inputVCFFile, final List<String> vcfSamples, Integer minVariantSamples, Integer gqThreshold, 
			final Double fractionSamplesPassing, List<String> ignoredChromosones, Log log, boolean filterToHetSNPs) {
		
		IOUtil.assertFileIsReadable(inputVCFFile);
		final VCFFileReader vcfReader = new VCFFileReader(inputVCFFile, false);
		log.info("Searching for variants with at least [" + minVariantSamples+ "] samples with the non-ref genotype");

		PeekableIterator<VariantContext> vcfIterator = SampleAssignmentVCFUtils.getVCFIterator(vcfReader, vcfSamples, false, gqThreshold, fractionSamplesPassing, ignoredChromosones, log);
		// if the genotype quality is disabled, be extra careful and filter out flip-snps (A/T, C/G SNPs that can be easily messed up in VCFs)
		// this is done by default in SampleAssignmentVCFUtils.getVCFIterator

		// filter monomorphic SNPs.
		vcfIterator = new PeekableIterator<>(new MonomorphicVariantContextFilter(vcfIterator, vcfSamples));
		// filter problematic chromosomes as defined by the user, usually the sex chromosomes.
		vcfIterator = new PeekableIterator<>(new CommonVariantContextFilter(vcfIterator, vcfSamples, minVariantSamples));
		// add a filter for SNPs that are singletons if you're in private SNP mode.
		if (filterToHetSNPs)
			vcfIterator = new PeekableIterator<>(new VariantContextSingletonFilter(vcfIterator, filterToHetSNPs));

		return vcfIterator;
	}
	
	/**
	 * Start with the list of all samples in the VCF.  If there's a file containing a subset of samples, 
	 * make sure all elements of the list are samples in the VCF.  If not, return null.
	 * @param vcfReader
	 * @param sampleFile
	 * @return The list of samples to analyze.  If only the vcfReader is supplied, the result is the samples in the VCF.  If the sampleFile is supplied, return that subset of samples if valid.  Return null if the list is invalid.
	 */
	public static List<String> getFinalSamplelist (VCFFileReader vcfReader, File sampleFile) {
		List<String> vcfSamples  = new ArrayList<>(vcfReader.getFileHeader().getSampleNameToOffset().keySet());
		// if a sample list is passed in, restrict list to those samples.
		if (sampleFile!=null)  {
			IOUtil.assertFileIsReadable(sampleFile);
			vcfSamples=ParseBarcodeFile.readCellBarcodeFile(sampleFile);
			List<String> validVcfSamples = SampleAssignmentVCFUtils.validateSampleNamesInVCF(vcfReader, vcfSamples, log);
			vcfSamples.removeAll(validVcfSamples);
			if (vcfSamples.size()>0) {
				log.info("Samples found in sample list but not VCF, quitting " + vcfSamples.toString()+"");
				return null;
			} else
				vcfSamples=validVcfSamples;
		}
		return (vcfSamples);
	}
	
	public static File getTempVCFFile (File tempDir) {
		File outTempVCF=null;		
		try {
			outTempVCF=File.createTempFile("tmp_vcf_", ".txt.gz", tempDir);
			outTempVCF.deleteOnExit();
			log.info("Writing temp VCF to " + outTempVCF.getAbsolutePath());
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		return (outTempVCF);

	}
	
	public static char getAltBase (final VariantContext vc) {
		Allele alt = vc.getAltAlleleWithHighestAlleleCount();

		char altBase='N';
		if (alt!=null) {
			byte [] altBases = alt.getBases();
			if (altBases.length>0)
				altBase=StringUtil.byteToChar(altBases[0]);
		}
		return altBase;
	}

	public static int compareRecords (final SNPGenomicBasePileUp p, final VariantContext vc, final SAMSequenceDictionary dict) {
		Interval sgpI = p.getSNPInterval();
		Interval sgpVC = new Interval (vc.getContig(), vc.getStart(), vc.getEnd());
		int cmp = IntervalTagComparator.compare(sgpI, sgpVC, dict);
		return cmp;
	}

	public static int compareRecords (final VariantContext vc1, final VariantContext vc2, final SAMSequenceDictionary dict) {
		if (vc1==null) return 1;
		if (vc2==null) return -1;
		Interval i1 = new Interval (vc1.getContig(), vc1.getStart(), vc1.getEnd());
		Interval i2 = new Interval (vc2.getContig(), vc2.getStart(), vc2.getEnd());
		int cmp = IntervalTagComparator.compare(i1, i2, dict);
		return cmp;
	}


}
