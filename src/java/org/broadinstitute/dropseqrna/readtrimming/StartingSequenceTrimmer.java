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
package org.broadinstitute.dropseqrna.readtrimming;

import htsjdk.samtools.util.SequenceUtil;
import htsjdk.samtools.util.StringUtil;
import picard.util.ClippingUtility;

/**
 * Abstract class that encapsulates the search for a starting sequence, and abstracts the algorithm for
 * determining how many mismatches is too many.  This class searches for the sequence starting at the end of the read,
 * and if the full sequence is found, but there are additional bases that precede the full sequence, that is still
 * considered a match.
 */
public abstract class StartingSequenceTrimmer {
 protected final int minMatch;
 protected final byte[] adapterSequence;

 public StartingSequenceTrimmer(final String sequence, final int minMatch) {
  this.minMatch = minMatch;
  this.adapterSequence = StringUtil.stringToBytes(sequence);
 }

 TrimResult trim(final String readBases) {
  // Inspired by Picard ClippingUtility.findIndexOfClipSequence
  if (readBases.length() < minMatch) {
   // If read is too short, can't possibly match.
   return StartingSequenceTrimmer.NO_TRIM;
  }
  final byte[] read = StringUtil.stringToBytes(readBases);
  // There may be multiple copies of the searched-for sequence.  We want to find the rightmost one, so start at the
  // end and walk backwards.
  int startIndex = -1;
  int endIndex;
  for (endIndex = read.length; endIndex >= minMatch; --endIndex) {
   startIndex = testMatch(read, endIndex);
   if (startIndex != -1) break;
  }

  if (startIndex == ClippingUtility.NO_MATCH) {
   return StartingSequenceTrimmer.NO_TRIM;
  }
  final boolean completelyTrimmed = endIndex == read.length;
  return new TrimResult(startIndex, endIndex, completelyTrimmed);

 }

 /**
  * Compare read to adapterSequence, walking backward, until read is exhausted, adapterSequence is exhausted, or
  * too many mismatches.
  * @param read bases to be tested
  * @param endIndex start test at endIndex-1
  * @return index in read where match starts, or -1 if no match.
  */
 private int testMatch(byte[] read, int endIndex) {
  final int possibleMatchLength = Math.min(endIndex, adapterSequence.length);
  final MismatchCounter mismatchCounter = getMismatchCounter();
  mismatchCounter.initialize(possibleMatchLength);
  int mismatches = 0;
  int readIndex = endIndex - 1;
  for (int adapterIndex = adapterSequence.length - 1;
       readIndex >= 0 && adapterIndex >= 0; --readIndex, --adapterIndex ) {
   if (!SequenceUtil.isNoCall(adapterSequence[adapterIndex]) &&
           !SequenceUtil.basesEqual(adapterSequence[adapterIndex], read[readIndex]) &&
           mismatchCounter.countMismatch()) {
    return -1;
   }
  }
  // This was decremented past the best position, so add 1
  return readIndex + 1;

 }
 protected abstract MismatchCounter getMismatchCounter();

 protected interface MismatchCounter {
  abstract void initialize(final int length);
  abstract boolean countMismatch();
 }

 public static class TrimResult {
  /**
   * 0-based position in read where starting sequence starts.  If starting sequence starts at beginning of the read,
   * this is 0.
   */
  public final int startPosition;

  /**
   * 0-based position in read after end of starting sequence.  E.g. if 10 bases of starting sequence are matched,
   * and start at the beginning of the read, then endPosition == 10.  If starting sequence is not found, endPosition==0.
   */
  public final int endPosition;

  /**
   * true if endPosition == length of read searched
   */
  public final boolean completelyTrimmed;

  public TrimResult(int startPosition, int endPosition, boolean completelyTrimmed) {
   this.startPosition = startPosition;
   this.endPosition = endPosition;
   this.completelyTrimmed = completelyTrimmed;
  }

  public boolean wasTrim() {
   return endPosition > 0;
  }
 }
 public static TrimResult NO_TRIM = new TrimResult(0, 0, false);
}
