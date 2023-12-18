/*
 * MIT License
 *
 * Copyright 2022 Broad Institute
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

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import org.broadinstitute.dropseqrna.utils.readpairs.ReadPair;

import java.util.Map;

/**
 * Iterate over queryname-sorted paired-read iterator and perform some validations.
 */
public class PairedSamRecordIterator
implements CloseableIterator<ReadPair> {
 private final CloseableIterator<SAMRecord> iterator;

 public PairedSamRecordIterator(CloseableIterator<SAMRecord> iterator) {
  this.iterator = iterator;
 }

 @Override
 public void close() {
  iterator.close();
 }

 @Override
 public boolean hasNext() {
  return iterator.hasNext();
 }

 @Override
 public ReadPair next() {
  SAMRecord r1 = iterator.next();
  // Checking for all sorts of conditions that will never happen.
  if (!r1.getReadPairedFlag()) {
   throw new SAMException("Unexpected unpaired read " + r1);
  }
  if (!iterator.hasNext()) {
   throw new SAMException("Premature EOF without mate for " + r1);
  }
  SAMRecord r2 = iterator.next();
  if (!r1.getReadName().equals(r2.getReadName())) {
   throw new SAMException(String.format("Mate not found for %s; Found %s instead of mate.", r1, r2));
  }
  if (!r2.getReadPairedFlag()) {
   throw new SAMException("Unexpected unpaired read " + r2);
  }
  if (r1.getFirstOfPairFlag()) {
   if (r2.getFirstOfPairFlag()) {
    throw new SAMException("Both reads are first of pair " + r1);
   }
  } else if (r2.getSecondOfPairFlag()) {
   throw new SAMException("Both reads are second of pair " + r1);
  }
  return new ReadPair(r1, r2);
 }
}
