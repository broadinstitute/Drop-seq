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
import java.util.Collection;
import java.util.List;
import java.util.Set;

import org.broadinstitute.dropseqrna.utils.IntervalTagComparator;
import org.broadinstitute.dropseqrna.vcftools.SampleAssignmentVCFUtils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;

public class CensusSeqUtils {


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

	//TODO: delete because unused?
	/*
	public static VariantContextWriter getVCFWriter (final File out, final VCFFileReader vcfReader, final VCFFileReader siteReader, final List<String> vcfSamples) {
		IOUtil.assertFileIsWritable(out);
		VariantContextWriter vcfWriter = SampleAssignmentVCFUtils.getVCFWriter(vcfReader, out);

		List<VCFHeader> headers = new ArrayList<>();
        headers.add(vcfReader.getFileHeader());
        headers.add(siteReader.getFileHeader());
        final VCFHeader mergedHeader = new VCFHeader(VCFUtils.smartMergeHeaders(headers, false), Collections.EMPTY_LIST);
        Set<VCFHeaderLine> mergedMetaData = mergedHeader.getMetaDataInInputOrder();
		VCFHeader newHeader = new VCFHeader(mergedMetaData, vcfSamples); // set up the new header with the restricted list of samples.
		vcfWriter.writeHeader(newHeader);
		return (vcfWriter);

	}
	*/

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
