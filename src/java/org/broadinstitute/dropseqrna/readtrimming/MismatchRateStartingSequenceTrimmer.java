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

/**
 * Starting sequence trimmer that has an error rate as opposed to a fixed number of errors.
 */
public class MismatchRateStartingSequenceTrimmer
        extends StartingSequenceTrimmer {
 private final double maxErrorRate;
 private final MismatchRateCounter mismatchRateCounter = new MismatchRateCounter();

 public MismatchRateStartingSequenceTrimmer(final String sequence, final int minMatch, final double maxErrorRate) {
  super(sequence, minMatch);
  this.maxErrorRate = maxErrorRate;
 }

 @Override
 protected MismatchCounter getMismatchCounter() {
  return mismatchRateCounter;
 }

 private class MismatchRateCounter
 implements MismatchCounter {
  private int mismatchesAllowed = -1;
  private int mismatches = -1;

  @Override
  public void initialize(int length) {
   mismatchesAllowed = (int)(length * maxErrorRate);
   mismatches = 0;
  }

  @Override
  public boolean countMismatch() {
   return ++mismatches > mismatchesAllowed;
  }
 }
}
