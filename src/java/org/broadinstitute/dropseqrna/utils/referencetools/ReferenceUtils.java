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
package org.broadinstitute.dropseqrna.utils.referencetools;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;

import java.util.ArrayList;
import java.util.List;

public class ReferenceUtils {

	public static String getSequence (byte [] fastaRefBases, Interval interval) {
		int startBase=interval.getStart();
		int endBase=interval.getEnd();
		byte [] bases=getSubArray (fastaRefBases, startBase-1, endBase-1);
		StringBuilder b= new StringBuilder();
		String baseString = new String(bases);
		b.append(baseString);
		return (b.toString());
	}
	
	private static byte [] getSubArray (byte [] input, int startPos, int endPos) {
		byte [] result = new byte[(endPos-startPos)+1];
		
		for (int i=0; i<=endPos-startPos; i++) {
			result[i]=input[i+startPos];
		}
		return result;
	}
	
	public static List<Interval> getIntervalsForChr (IntervalList intervals, String chromosome) {
		List<Interval> result = new ArrayList<Interval>();
		for (Interval i: intervals) {
			if (i.getContig().equals(chromosome)) {
				result.add(i);
			}
		}
		return (result);	
	}
	
	
}
