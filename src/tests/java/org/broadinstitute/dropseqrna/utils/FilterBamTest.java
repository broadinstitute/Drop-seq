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
package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.*;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.samtools.util.zip.DeflaterFactory;
import org.testng.Assert;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import picard.sam.ValidateSamFile;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

// Note that there is some weirdness with IntelDeflater on Mac OS that sometimes causes SEGV.  This test therefore
// tried to disable IntelDeflater as much as possible when running on Mac.  The problem does not always happen if just
// a single test is run.  Perhaps it has to do with IntelDeflater being used multiple times in same process.
public class FilterBamTest {

	// TO GENERATE COUNTS OF CONTIGS FOR HUMAN/MOUSE...
	// samtools view human_mouse_smaller.bam |cut -f 3 |sort |uniq -c > human_mouse_smaller.contig_counts.txt
	private static File TEST_DATA_DIR = new File("testdata/org/broadinstitute/dropseq/utils");
	private static File ORGANISM_INPUT_FILE=new File (TEST_DATA_DIR,"human_mouse_smaller.bam");
	private static File BAM_LIST = new File(TEST_DATA_DIR, "human_mouse_x.bam_list");
	private static List<File> BAM_FILES = Arrays.asList(
			new File(TEST_DATA_DIR, "human_mouse_x.1.sam"),
			new File(TEST_DATA_DIR, "human_mouse_x.2.sam")
	);

	private static final long NUM_HUMAN_READS = 136859;
	public static final long NUM_MOUSE_READS = 76036;
	public static final long TOTAL_READS = 248661;

	@Test
	public void testOrganismFiltering () throws IOException {
		FilterBam f = new FilterBam();
		f.INPUT= Collections.singletonList(ORGANISM_INPUT_FILE);
		f.SUMMARY=File.createTempFile("paired_input", ".summary.txt");
		f.OUTPUT=File.createTempFile("paired_input", ".bam");
		f.OUTPUT.deleteOnExit();
		f.SUMMARY.deleteOnExit();
		f.REF_SOFT_MATCHED_RETAINED=Arrays.asList("HUMAN");
		// need to set this to an empty list, as that's what the constructor for FilterBam would do.
		f.STRIP_REF_PREFIX=new ArrayList<>();
		f.doWork();
		
		List<FilteredReadsMetric> metrics = MetricsFile.readBeans(f.SUMMARY);
		Assert.assertEquals(NUM_HUMAN_READS, metrics.get(0).READS_ACCEPTED);
		Assert.assertEquals(TOTAL_READS - NUM_HUMAN_READS, metrics.get(0).READS_REJECTED);
		
		f.REF_SOFT_MATCHED_RETAINED=Arrays.asList("MOUSE");
		Assert.assertEquals(f.doWork(), 0);
		metrics = MetricsFile.readBeans(f.SUMMARY);
		Assert.assertEquals(NUM_MOUSE_READS, metrics.get(0).READS_ACCEPTED);
		Assert.assertEquals(TOTAL_READS - NUM_MOUSE_READS, metrics.get(0).READS_REJECTED);
		
		
	}

	@Test(dataProvider = "thresholdTestDataProvider")
	public void testPassingThreshold(final boolean fractionalThreshold) throws IOException  {
		filterWithThreshold(fractionalThreshold? (NUM_HUMAN_READS - 1)/(double)TOTAL_READS: NUM_HUMAN_READS);
	}

	@Test(dataProvider = "thresholdTestDataProvider", expectedExceptions = {RuntimeException.class})
	public void testFailingThreshold(final boolean fractionalThreshold) throws IOException  {
		filterWithThreshold(fractionalThreshold? (NUM_HUMAN_READS + 1)/(double)TOTAL_READS: NUM_HUMAN_READS + 1);
	}
	private void filterWithThreshold(final double threshold) throws IOException {
		FilterBam f = new FilterBam();
		f.INPUT=Collections.singletonList(ORGANISM_INPUT_FILE);;
		f.OUTPUT=File.createTempFile("paired_input", ".bam");
		f.OUTPUT.deleteOnExit();
		f.REF_SOFT_MATCHED_RETAINED=Arrays.asList("HUMAN");
		// need to set this to an empty list, as that's what the constructor for FilterBam would do.
		f.STRIP_REF_PREFIX=new ArrayList<>();
		f.PASSING_READ_THRESHOLD = threshold;
		Assert.assertEquals(f.doWork(), 0);
	}

	@DataProvider(name="thresholdTestDataProvider")
	public Object[][] thresholdTestDataProvider() {
		return new Object[][] {
				{true},
				{false}
		};
	}
	
	@Test(enabled=true, groups = { "dropseq","transcriptome" })
	public void testCigarFilter() {
		FilterBam b = new FilterBam();
		b.SUM_MATCHING_BASES=40;


		SAMRecord r = new SAMRecord(null);
		r.setCigarString("50M");
		boolean rejFlag = b.rejectOnCigar(r);
		Assert.assertFalse(rejFlag);

		boolean result = b.filterRead(r);
		Assert.assertFalse(result);


		r.setCigarString("20M");
		rejFlag = b.rejectOnCigar(r);
		Assert.assertTrue(rejFlag);

		result = b.filterRead(r);
		Assert.assertTrue(result);

		r.setCigarString("20S25M30H");
		rejFlag = b.rejectOnCigar(r);
		Assert.assertTrue(rejFlag);

		result = b.filterRead(r);
		Assert.assertTrue(result);

		r.setCigarString("20M10S20M");
		rejFlag = b.rejectOnCigar(r);
		Assert.assertFalse(rejFlag);

		result = b.filterRead(r);
		Assert.assertFalse(result);




	}

	@Test(enabled=true, groups = { "dropseq","transcriptome" })
	public void testNonPrimary() {
		FilterBam b = new FilterBam();
		b.RETAIN_ONLY_PRIMARY_READS=true;

		SAMRecord r = new SAMRecord(null);
		r.setSecondaryAlignment(true);
		boolean result = b.filterRead(r);
		Assert.assertTrue(result);

		r = new SAMRecord(null);
		r.setSecondaryAlignment(false);
		result = b.rejectNonPrimaryReads(r);
		Assert.assertFalse(result);

		b.RETAIN_ONLY_PRIMARY_READS=false;
		r.setSecondaryAlignment(true);
		result = b.filterRead(r);
		Assert.assertFalse(result);
	}

	@Test(enabled=true, groups = { "dropseq", "transcriptome" })
	public void testMapQuality () {
		FilterBam b = new FilterBam();
		b.MINIMUM_MAPPING_QUALITY=10;
		SAMRecord r = new SAMRecord(null);
		r.setMappingQuality(12);
		Assert.assertFalse(b.filterRead(r));

		r.setMappingQuality(9);
		Assert.assertTrue(b.filterRead(r));

		b.MINIMUM_MAPPING_QUALITY=null;
		Assert.assertFalse(b.filterRead(r));
	}

	@Test(enabled=true, groups = { "dropseq", "transcriptome" })
	public void testSoftMatch () {
		FilterBam b = new FilterBam();
		List<String> retained = new ArrayList<>();
		retained.add("HUMAN");
		b.REF_SOFT_MATCHED_RETAINED=retained;
		b.buildPatterns();
		SAMRecord r = new SAMRecord(null);
		r.setReferenceName("HUMAN_CHR1");
		Assert.assertFalse(b.filterRead(r));
		r.setReferenceName("CHIMP_CHR1");
		Assert.assertTrue(b.filterRead(r));

		retained.add("CHIMP");
		b.REF_SOFT_MATCHED_RETAINED=retained;
		b.buildPatterns();
		Assert.assertFalse(b.filterRead(r));

		r.setReferenceName("HUMAN_CHR1");
		List<String> matches = Arrays.asList("HUMAN");
		Assert.assertTrue(b.softMatchReference(matches, r));
		matches = Arrays.asList("MOUSE");
		Assert.assertFalse(b.softMatchReference(matches, r));


	}

	@Test(enabled=true, groups = { "dropseq", "transcriptome" })
	public void testExactMatch () {
		FilterBam b = new FilterBam();
		List<String> retained = new ArrayList<>();
		retained.add("HUMAN_CHR1");
		b.REF_HARD_MATCHED_RETAINED=retained;
		b.buildPatterns();
		SAMRecord r = new SAMRecord(null);
		r.setReferenceName("HUMAN_CHR1");
		Assert.assertFalse(b.filterRead(r));
		Assert.assertTrue(b.exactMatchReference(retained, r));
		r.setReferenceName("CHIMP_CHR1");
		Assert.assertTrue(b.filterRead(r));
		Assert.assertFalse(b.exactMatchReference(retained, r));




	}

	@Test(enabled=true, groups = { "dropseq", "transcriptome" })
	public void testRejectOnTags () {
		FilterBam b = new FilterBam();
		b.TAG_REJECT_COMBINE_FLAG="UNION";
		List<String> tags = new ArrayList<>();
		tags.add("XE");
		tags.add("XZ");
		b.TAG_REJECT=tags;
		SAMRecord r = new SAMRecord(null);
		r.setAttribute("XE", "Foo");
		boolean result = b.filterRead(r);
		Assert.assertTrue(result);

		r = new SAMRecord(null);
		r.setAttribute("ZZ", "Foo");
		result = b.filterRead(r);
		Assert.assertFalse(result);

		// if you need to see 2 tags and only see 1, then you don't reject the read.
		b.TAG_REJECT_COMBINE_FLAG="INTERSECT";
		result = b.rejectOnTags(tags, r);
		Assert.assertFalse(result);
	}

	@Test(enabled=true, groups = { "dropseq", "transcriptome" })
	public void testAcceptOnTags () {
		FilterBam b = new FilterBam();
		List<String> tags = new ArrayList<>();
		tags.add("XE");
		tags.add("XZ");
		b.TAG_RETAIN=tags;
		SAMRecord r = new SAMRecord(null);
		r.setAttribute("XE", "Foo");
		boolean result = b.filterRead(r);
		Assert.assertFalse(result);

		r = new SAMRecord(null);
		r.setAttribute("QQ", "Foo");
		result = b.filterRead(r);
		Assert.assertTrue(result);
	}


	private SAMRecord makeAlignment(final SAMFileHeader header,
									final int referenceIndex,
									final int alignmentStart,
									final int mateReferenceIndex,
									final int mateAlignmentStart,
									final String readName,
									final String readGroupId) {

		final String readBases = "ACGT";
		final String baseQualities = "2222";
		final SAMRecord rec = new SAMRecord(header);
		rec.setReferenceIndex(referenceIndex);
		rec.setAlignmentStart(alignmentStart);
		rec.setReadString(readBases);
		rec.setBaseQualityString(baseQualities);
		rec.setCigarString(readBases.length() + "M");
		rec.setReadName(readName);
		rec.setAttribute("RG", readGroupId);
		rec.setAttribute("NM", 0);
		rec.setReadPairedFlag(true);
		rec.setFirstOfPairFlag(true);
		if (mateReferenceIndex == SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX)
			rec.setMateUnmappedFlag(true);
		else {
		    rec.setMateReferenceIndex(mateReferenceIndex);
		    rec.setMateAlignmentStart(mateAlignmentStart);
        }
		return rec;
	}

	private File createInputSam() throws IOException {
		final File inputSamFile = File.createTempFile("FilterBamTest.", ".sam");
		final File referenceFasta = File.createTempFile("FilterBamTest.", ".fasta");
		inputSamFile.deleteOnExit();
		referenceFasta.deleteOnExit();

		final SAMSequenceDictionary sd = new SAMSequenceDictionary();
		sd.addSequence(new SAMSequenceRecord("MOUSE_m1", 1000000));
		sd.addSequence(new SAMSequenceRecord("HUMAN_h1", 1000000));
		sd.addSequence(new SAMSequenceRecord("MOUSE_m2", 1000000));
		sd.addSequence(new SAMSequenceRecord("HUMAN_h2", 1000000));
		final SAMFileHeader header = new SAMFileHeader(sd);

        SAMReadGroupRecord readGroupRecord = new SAMReadGroupRecord("1");
        readGroupRecord.setSample("FREE_SAMPLE");
        readGroupRecord.setPlatform("ILLUMINA");
        List<SAMReadGroupRecord> readGroups = new ArrayList();
        readGroups.add(readGroupRecord);
        header.setReadGroups(readGroups);

        final SAMFileWriter writer = new SAMFileWriterFactory().makeWriter(header, false, inputSamFile, referenceFasta);
		for (int i = 0; i < sd.getSequences().size(); ++i)
			writer.addAlignment(makeAlignment(header, i, i+1, SAMRecord.NO_ALIGNMENT_REFERENCE_INDEX, 0,
                    "read" + i, readGroupRecord.getId()));

		// Add some read with mate info, some that have same-species mate, and some that have other-species mate
        writer.addAlignment(makeAlignment(header, 0, 1, 2, 1,
                "mouse_mouse", readGroupRecord.getId()));
        writer.addAlignment(makeAlignment(header, 1, 1, 3, 1,
                "human_human", readGroupRecord.getId()));
        writer.addAlignment(makeAlignment(header, 0, 1, 1, 1,
                "mouse_human", readGroupRecord.getId()));
        writer.addAlignment(makeAlignment(header, 1, 1, 2, 1,
                "human_house", readGroupRecord.getId()));

		writer.close();
		Assert.assertTrue(validateSamFile(inputSamFile));
		return inputSamFile;
	}

	private boolean validateSamFile(final File bamFile) {
		final int ret = new ValidateSamFile().instanceMain(new String[]{
				"INPUT=" + bamFile.getAbsolutePath(),
				"OUTPUT=/dev/stdout",
                // reads with mate info are created, but we don't bother to create the mates, so ignore this error.
                "IGNORE=MATE_NOT_FOUND",
				"USE_JDK_DEFLATER=" + TestUtils.isMacOs()});
		return ret == 0;
	}

	@Test(dataProvider = "testEndToEndDataProvider")
	public void testEndToEnd(final String organism1, final String organism2, final boolean stripPrefix,
                             final boolean dropSequences, final boolean soft, final boolean retain) throws IOException {
		final File inputSam = createInputSam();
		final File outputSam = File.createTempFile("FilterBamTest." +
                String.format("stripPrefix_%s;dropSequences_%s;soft_%s;retain_%s", stripPrefix, dropSequences, soft, retain) + ".",
                ".sam");
		outputSam.deleteOnExit();
		final ArrayList<String> args = new ArrayList<>(Arrays.asList(
		        "INPUT=" + inputSam.getAbsolutePath(),
                "OUTPUT=" + outputSam.getAbsolutePath(),
                "DROP_REJECTED_REF=" + dropSequences
				));
		if (TestUtils.isMacOs()) {
			args.add("USE_JDK_DEFLATER=true");
		}
		if (stripPrefix) {
		    args.add("STRIP_REF_PREFIX=" + organism1);
		    args.add("STRIP_REF_PREFIX=" + organism2);
        }
        final String ref_filter;
        if (soft) {
            if (retain)
				ref_filter = "REF_SOFT_MATCHED_RETAINED";
			else
				ref_filter = "REF_SOFT_MATCHED_REJECTED";
        } else if (retain)
			ref_filter = "REF_HARD_MATCHED_RETAINED";
		else
			ref_filter = "REF_HARD_MATCHED_REJECTED";
        args.add(ref_filter + "=" + organism1);
        final int ret = new FilterBam().instanceMain(args.toArray(new String[args.size()]));
        Assert.assertEquals(0, ret);
        Assert.assertTrue(validateSamFile(outputSam));
	}

	@DataProvider(name="testEndToEndDataProvider")
    public Object[][] testEndToEndDataProvider() {
	    final ArrayList<Object[]> ret = new ArrayList<>();
	    for (final boolean stripPrefix : Arrays.asList(true, false))
			for (final boolean dropSequences: Arrays.asList(true, false))
				for (final boolean soft: Arrays.asList(true, false))
					for (final boolean retain : Arrays.asList(true, false))
						ret.add(new Object[]{"HUMAN_", "MOUSE_", stripPrefix, dropSequences, soft, retain});
        return ret.toArray(new Object[ret.size()][]);
    }

	@Test(dataProvider = "testMultipleBamInputDataProvider")
	public void testMultipleBamInput(final List<File> inputBamsOrBamList) throws IOException {
		FilterBam f = new FilterBam();
		f.INPUT= inputBamsOrBamList;
		f.SUMMARY=TestUtils.getTempReportFile("human_mouse.x", ".summary.txt");
		f.OUTPUT=TestUtils.getTempReportFile("paired_input", ".bam");
		f.REF_SOFT_MATCHED_RETAINED=Arrays.asList("HUMAN");
		// need to set this to an empty list, as that's what the constructor for FilterBam would do.
		f.STRIP_REF_PREFIX=new ArrayList<>();
		Assert.assertEquals(f.doWork(), 0);
		// Just count the reads in the input files because we know no reads are filtered out
		long numReads = BAM_FILES.stream().mapToLong(TestUtils::countSamRecords).sum();
		List<FilteredReadsMetric> metrics = MetricsFile.readBeans(f.SUMMARY);
		Assert.assertEquals(numReads, metrics.get(0).READS_ACCEPTED);
	}
	@DataProvider(name="testMultipleBamInputDataProvider")
	public Object[][] testMultipleBamInputDataProvider() {
		return new Object[][] {
				{BAM_FILES},
				{Collections.singletonList(BAM_LIST)}
		};
	}

	private DeflaterFactory deflaterFactory;
    @BeforeMethod
	public void setDeflaterFactory() {
    	if (TestUtils.isMacOs()) {
			deflaterFactory = BlockCompressedOutputStream.getDefaultDeflaterFactory();
			BlockCompressedOutputStream.setDefaultDeflaterFactory(new DeflaterFactory());
		}
	}
	@AfterMethod
	public void restoreDeflaterFactory() {
		if (TestUtils.isMacOs()) {
			if (deflaterFactory != null) {
				BlockCompressedOutputStream.setDefaultDeflaterFactory(deflaterFactory);
			}
		}
	}
}
