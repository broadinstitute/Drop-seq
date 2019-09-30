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

package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.CloserUtil;
import org.testng.Assert;
import picard.util.TabbedInputParser;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;

/**
 * Helper methods for
 * @author nemesh
 *
 */
public class TestUtils {

	/**
	 * Test if two text files are the same, ignoring "#" character lines.
	 * @param expected
	 * @param actual
	 * @return
	 */
	public static boolean testFilesSame (final File expected, final File actual) {
		TabbedInputParser e = new TabbedInputParser(true, expected);
		TabbedInputParser a = new TabbedInputParser(true, actual);
		while (e.hasNext() && a.hasNext()) {
			e.next();
			a.next();
			String le = e.getCurrentLine();
			String la = a.getCurrentLine();
			if (!le.equals(la)) {
				CloserUtil.close(e);
				CloserUtil.close(a);
				return false;
			}
		}
		CloserUtil.close(e);
		CloserUtil.close(a);
		return true;
	}

	// Copied roughly from Picard CompareMetrics
	public static boolean testMetricsFilesEqual(final File expected, final File actual) throws FileNotFoundException {
		final MetricsFile<?, ?> metricsA = new MetricsFile();
		final MetricsFile<?, ?> metricsB = new MetricsFile();
		metricsA.read(new FileReader(expected));
		metricsB.read(new FileReader(actual));
		return metricsA.areMetricsEqual(metricsB) && metricsA.areHistogramsEqual(metricsB);
	}

	public static void assertSamFilesSame(final File actual, final File expected) {
		final SamReader expectedReader = SamReaderFactory.makeDefault().open(expected);
		final SamReader actualReader = SamReaderFactory.makeDefault().open(actual);
		Assert.assertEquals(expectedReader.getFileHeader(), actualReader.getFileHeader());
		final SAMRecordIterator expectedIterator = expectedReader.iterator();
		final SAMRecordIterator actualIterator = actualReader.iterator();
		while (expectedIterator.hasNext()) {
			if (!actualIterator.hasNext()) {
				Assert.fail(String.format("%s has more records than %s", expected, actual));
			} else {
				Assert.assertEquals(actualIterator.next(), expectedIterator.next());
			}
		}
		if (actualIterator.hasNext()) {
			Assert.fail(String.format("%s has fewer records than %s", expected, actual));
		}
	}

	public static boolean isMacOs() {
		return System.getProperty("os.name").toLowerCase().contains("mac");
	}
}
