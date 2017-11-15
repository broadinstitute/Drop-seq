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
package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

import java.io.File;
import java.util.Collections;

import junit.framework.Assert;

import org.testng.annotations.Test;

public class MapQualityProcessorTest {



	@Test(enabled = true)
	public void processReadTest() {
		String unmappedReadName="NS500217:67:H14GMBGXX:3:23408:5941:1275";
		String mappedReadName="NS500217:67:H14GMBGXX:1:22207:3769:12483";
		SAMRecord unmapped = getRecordFromBAM(unmappedReadName);
		SAMRecord mapped = getRecordFromBAM(mappedReadName);

		MapQualityFilteredIterator r1 = new MapQualityFilteredIterator(Collections.singletonList(unmapped).iterator(), 10, true);
		Assert.assertFalse(r1.hasNext());

        MapQualityFilteredIterator r2 = new MapQualityFilteredIterator(Collections.singletonList(mapped).iterator(), 10, true);
		Assert.assertTrue(r2.hasNext());
		Assert.assertEquals(mapped, r2.next());
        Assert.assertFalse(r2.hasNext());
	}


	// A slow way to get reads for testing.
	private SAMRecord getRecordFromBAM (final String readName) {
		File inFile = new File("testdata/org/broadinstitute/transcriptome/barnyard/5cell3gene.bam");
		SamReader reader = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE).open(inFile);
		for (SAMRecord r: reader)
			if (r.getReadName().equals(readName)) return (r);
		throw new IllegalArgumentException("Asked for a read not in the BAM!");
	}
}
