/*
 * MIT License
 *
 * Copyright 2021 Broad Institute
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

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;

public class ConvertTagToReadGroupTest {
	private static final File TESTDATA_DIR = new File("testdata/org/broadinstitute/dropseq/utils");
	private static final File ALIGNED_UNPAIRED_BAM = new File(TESTDATA_DIR, "N701_very_small.sam");
	private static final File NO_RG_BAM = new File(TESTDATA_DIR, "SequenceDictionaryIntersectionTest/no_chr.sam");
	private static final String TAG = "XC";
	private static final String BAD_TAG = "XX";
	private static final String SAMPLE = "s1";
	private static final String LIBRARY = "l1";
	private static final String PLATFORM_UNIT = "PU1";
	private static final String PLATFORM = "P1";

	@Test(dataProvider = "testBasicDataProvider")
	public void testBasic(final int READ_MQ, final boolean use_CELL_BC_FILE, final boolean optionalRgFields, final boolean allBarcodes) {
		final ConvertTagToReadGroup clp = new ConvertTagToReadGroup();
		clp.INPUT = ALIGNED_UNPAIRED_BAM;
		clp.OUTPUT = TestUtils.getTempReportFile("ConvertTagToReadGroupTest.", ".sam");
		clp.READ_MQ = READ_MQ;
		clp.SAMPLE_NAME = SAMPLE;
		final String library;
		final String platform_unit;
		final String platform;
		if (optionalRgFields) {
			library = LIBRARY;
			platform_unit = PLATFORM_UNIT;
			platform = PLATFORM;
			clp.LIBRARY_NAME = library;
			clp.PLATFORM_UNIT = platform_unit;
			clp.PLATFORM = platform;
		} else {
			SamReader in = SamReaderFactory.makeDefault().open(clp.INPUT);
			final SAMReadGroupRecord rg = in.getFileHeader().getReadGroups().get(0);
			CloserUtil.close(in);
			library = rg.getLibrary();
			platform_unit = rg.getPlatformUnit();
			platform = rg.getPlatform();
		}
		final ObjectCounter<String> barcodeCounts = DownsampleBamByTagTest.getTagCounts(Collections.singletonList(clp.INPUT), TAG, READ_MQ, false);
		final int numBarcodes = (allBarcodes ? barcodeCounts.getSize() : barcodeCounts.getSize() / 2);
		final Collection<String> expectedBarcodes;
		if (use_CELL_BC_FILE) {
			final List<String> originalBarcodes = new ArrayList<>(barcodeCounts.getKeys());
			Collections.shuffle(originalBarcodes);
			expectedBarcodes = originalBarcodes.subList(0, numBarcodes);
			clp.CELL_BC_FILE = TestUtils.getTempReportFile("ConvertTagToReadGroupTest.", ".cell_barcodes.txt");
			final PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(clp.CELL_BC_FILE));
			for (final String barcode : expectedBarcodes) {
				out.println(barcode);
			}
			out.close();
		} else {
			final List<String> originalBarcodes = barcodeCounts.getKeysOrderedByCount(true);
			expectedBarcodes = originalBarcodes.subList(0, numBarcodes);
			clp.NUM_CORE_BARCODES = numBarcodes;
		}
		Assert.assertEquals(clp.doWork(), 0);
		SamReader in = SamReaderFactory.makeDefault().open(clp.OUTPUT);
		final List<SAMReadGroupRecord> readGroups = in.getFileHeader().getReadGroups();
		Assert.assertEquals(readGroups.size(), expectedBarcodes.size());
		final Set<String> expectedRgIdSet = new HashSet<>(expectedBarcodes);
		final Set<String> expectedSampleSet = expectedBarcodes.stream().map(bc -> SAMPLE + ":" + bc).collect(Collectors.toSet());
		for (SAMReadGroupRecord rg : readGroups) {
			Assert.assertTrue(expectedRgIdSet.contains(rg.getId()));
			Assert.assertTrue(expectedSampleSet.contains(rg.getSample()));
			Assert.assertEquals(rg.getLibrary(), library);
			Assert.assertEquals(rg.getPlatformUnit(), platform_unit);
			Assert.assertEquals(rg.getPlatform(), platform);
		}
		for (final SAMRecord rec : in) {
			final String rgId = rec.getStringAttribute(SAMTag.RG.name());
			final String barcode = rec.getStringAttribute(TAG);
			Assert.assertEquals(rgId, barcode);
			Assert.assertTrue(expectedRgIdSet.contains(rgId));
		}
		CloserUtil.close(in);
	}

	@DataProvider(name = "testBasicDataProvider")
	public Object[][] testBasicDataProvider() {
		final boolean[] tf = { true, false };
		final int[] readMQs = { 20, 0 };
		final ArrayList<Object[]> ret = new ArrayList<>();
		for (final int READ_MQ : readMQs) {
			for (final boolean use_CELL_BC_FILE : tf) {
				for (final boolean optionalRgFields : tf) {
					for (final boolean allBarcodes : tf) {
						Object[] parameters = { READ_MQ, use_CELL_BC_FILE, optionalRgFields, allBarcodes };
						ret.add(parameters);
					}
				}
			}
		}
		return ret.toArray(new Object[0][]);
	}

	@Test
	public void testNoTag() {
		final ConvertTagToReadGroup clp = new ConvertTagToReadGroup();
		clp.INPUT = ALIGNED_UNPAIRED_BAM;
		clp.OUTPUT = TestUtils.getTempReportFile("ConvertTagToReadGroupTest.", ".sam");
		clp.READ_MQ = 0;
		clp.SAMPLE_NAME = SAMPLE;
		clp.NUM_CORE_BARCODES = 10;
		clp.CELL_BARCODE_TAG = BAD_TAG;
		Assert.assertEquals(clp.doWork(), 0);
		// Confirm that 0 reads are written
		SamReader in = SamReaderFactory.makeDefault().open(clp.OUTPUT);
		Assert.assertFalse(in.iterator().hasNext());
		CloserUtil.close(in);
	}

	@Test(expectedExceptions = { IllegalArgumentException.class })
	public void testNoRg() {
		final ConvertTagToReadGroup clp = new ConvertTagToReadGroup();
		clp.INPUT = NO_RG_BAM;
		clp.OUTPUT = TestUtils.getTempReportFile("ConvertTagToReadGroupTest.", ".sam");
		clp.READ_MQ = 0;
		clp.SAMPLE_NAME = SAMPLE;
		clp.NUM_CORE_BARCODES = 10;
		clp.CELL_BARCODE_TAG = BAD_TAG;
		Assert.assertNotEquals(clp.doWork(), 0);
	}
}
