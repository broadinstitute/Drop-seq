/*
 * MIT License
 *
 * Copyright 2023 Broad Institute
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

import htsjdk.samtools.*;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class UnsortedMergingSamRecordIteratorTest {

 // a couple of random SAM files
 private static final File[] INPUTS = {
         new File("testdata/org/broadinstitute/dropseq/barnyard/MarkChimericReads.input2.sam"),
         new File("testdata/org/broadinstitute/dropseq/readtrimming/N701.subset.tagged_filtered.sam")
 };

 @Test(dataProvider = "testBasicDataProvider")
 public void testBasic(boolean mergeHeaders) {
  final SamReaderFactory readerFactory = SamReaderFactory.makeDefault();
  final List<SamReader> samReaders = Arrays.stream(INPUTS).map(f -> readerFactory.open(f)).toList();
  final List<String> expectedReadNames = samReaders.stream().map(SamReader::iterator).
          flatMap(it -> it.stream()).map(SAMRecord::getReadName).toList();
  final SAMFileHeader header;
  if (mergeHeaders) {
   header = (new SamFileHeaderMerger(SAMFileHeader.SortOrder.unsorted, samReaders.stream().map(SamReader::getFileHeader).toList(),
           true)).getMergedHeader();
  } else {
   header = null;
  }
  final UnsortedMergingSamRecordIterator it = new UnsortedMergingSamRecordIterator(header,
          Arrays.stream(INPUTS).map(f -> readerFactory.open(f)).toList());
  final List<String> actualReadNames = it.stream().map(SAMRecord::getReadName).toList();
  Assert.assertEquals(actualReadNames, expectedReadNames);
 }

 @DataProvider(name = "testBasicDataProvider")
 private Object[][] testBasicDataProvider() {
  return new Object[][]{
          {true},
          {false}
  };
 }
}
