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
import java.util.Set;
import java.util.TreeSet;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.priv.barnyard.digitalallelecounts.sampleassignment.genomicpool.CensusSeqUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        summary = "Reads a VCF/VCF.gz/and extracts SNPs that are heterozygous for all listed samples",
        oneLineSummary = "Creates an interval file of variants from a VCF",
        programGroup = DropSeq.class)
/**
 * This is currently a very single purpose tool to set up GatherDigitalAlleleCounts by finding heterozygous intervals for one or more samples.
 * @author nemesh
 *
 */
public class CreateSnpIntervalFromVcf extends CommandLineProgram {

	private static final Log log = Log.getInstance(CreateSnpIntervalFromVcf.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input VCF.")
	public File INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The intervals file containing variant sites from the VCF that pass filters")
	public File OUTPUT;

	@Argument(doc = "The samples to that must be heterozygous for a variant to be included in the output.", optional=true)
	public Set<String> SAMPLE = new TreeSet<String>();

	@Argument(doc= "The minimum genotype quality for a variant across all samples.")
	public Integer GQ_THRESHOLD=0;

	@Argument(shortName = StandardOptionDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, doc = "The reference sequence dictionary.  This will be the header for the intervals output.", optional = true)
	public File SEQUENCE_DICTIONARY;

	@Argument(doc="Only output heterozygous SNPs for these samples.")
	public boolean HET_SNPS_ONLY=false;

	@Override
	protected int doWork() {
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsWritable(OUTPUT);
		IntervalList result = processData(this.INPUT, this.SEQUENCE_DICTIONARY, SAMPLE, this.GQ_THRESHOLD, this.HET_SNPS_ONLY);
		result.write(this.OUTPUT);
		return 0;
	}


	public IntervalList processData(final File vcfFile, final File sdFile, final Set<String> sample, int GQThreshold, final boolean hetSNPsOnly) {

		final VCFFileReader reader = new VCFFileReader(vcfFile, false);
		if (!CensusSeqUtils.GQInHeader(reader)) {
			GQThreshold=-1;
			log.info("Genotype Quality [GQ] not found in header.  Disabling GQ_THRESHOLD parameter");
		}
		
		final VCFHeader inputVcfHeader = new VCFHeader(reader.getFileHeader().getMetaDataInInputOrder());
		SAMSequenceDictionary sequenceDictionary = inputVcfHeader.getSequenceDictionary();
		Set<String> sampleListFinal = sample;
		if (sample==null || sample.isEmpty()) {
			ArrayList<String> s = reader.getFileHeader().getSampleNamesInOrder();
			sampleListFinal=new TreeSet<String>(s);
		}

		if (sdFile != null)
			sequenceDictionary = getSequenceDictionary(sdFile);

		final ProgressLogger progress = new ProgressLogger(this.log, 500000);

		final SAMFileHeader samHeader = new SAMFileHeader();
		samHeader.setSequenceDictionary(sequenceDictionary);
		IntervalList result = new IntervalList(samHeader);

		// Go through the input, find sites we want to keep.
		final PeekableIterator<VariantContext> iterator = new PeekableIterator<>(reader.iterator());

		validateRequestedSamples (iterator, sampleListFinal);

		while (iterator.hasNext()) {
			final VariantContext site = iterator.next();
			progress.record(site.getContig(), site.getStart());
			// for now drop any filtered site.
			if (site.isFiltered())
				continue;
			// move onto the next record if the site is not a SNP or the samples aren't all heterozygous.
			if (!site.isSNP())
				continue;
			if (!sitePassesFilters(site, sampleListFinal, GQThreshold, hetSNPsOnly))
				continue;
			Interval varInt = new Interval(site.getContig(), site.getStart(),
					site.getEnd(), true, site.getID());

			// final Interval site = findHeterozygousSites(full, SAMPLE);
			result.add(varInt);

		}

		CloserUtil.close(iterator);
		CloserUtil.close(reader);
		return (result);
	}

	private void validateRequestedSamples (final PeekableIterator<VariantContext> iterator, final Set<String> sample) {
		if (iterator.hasNext()) {
			final VariantContext site = iterator.peek();
			Set<String> vCFSamples = site.getSampleNames();
			for (String s: sample)
				if (!vCFSamples.contains(s))
					throw new IllegalArgumentException("VCF doesn't have the requested sample " + s);
		}
	}

	private static SAMSequenceDictionary getSequenceDictionary(final File sd) {
		SamReader reader = SamReaderFactory.makeDefault().open(sd);
		SAMSequenceDictionary dict = reader.getFileHeader().getSequenceDictionary();
		// final SAMFileReader r = new SAMFileReader(sd);
		// SAMSequenceDictionary dict = r.getFileHeader().getSequenceDictionary();
		CloserUtil.close(reader);

		return dict;
	}

	/** Are all the samples listed heterozygous for this variant, and pass GQ threshold?*/
	private boolean sitePassesFilters(final VariantContext ctx, final Set<String> samples, final int GQThreshold, final boolean hetSNPsOnly) {

		boolean flag = true;
		for (String sample : samples) {
			Genotype g = ctx.getGenotype(sample);
			if (g.getGQ()<GQThreshold || (!g.isHet() && hetSNPsOnly))
			 return (false);
		}
		return (flag);
	}

	// Stock main method
	public static void main(final String[] args) {
		new CreateSnpIntervalFromVcf().instanceMainWithExit(args);
	}

}
