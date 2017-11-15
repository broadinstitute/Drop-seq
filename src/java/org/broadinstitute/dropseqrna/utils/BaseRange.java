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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class BaseRange {

	private int start;
	private int end;

	public BaseRange(final int start, final int end) {
		this.start= start;
		this.end=end;
	}

	public int getSize () {
		int r = this.getEnd()-this.getStart()+1;
		return (r);
	}

	public static int getTotalRangeSize (final String baseRange) {
		List<BaseRange> result = parseBaseRange(baseRange);
		int totalSize=0;
		for (BaseRange b: result) {
			int s = b.getSize();
			totalSize+=s;
		}
		return totalSize;
	}

	public static List<BaseRange> parseBaseRange(final String baseRange) {
		String [] split = baseRange.split(":");
		List<BaseRange> result = new ArrayList<BaseRange>(split.length);
		for (String s: split) {
			BaseRange r = parseSingleBaseRange(s);
			result.add(r);
		}
		return result;
	}

	public static BaseRange parseSingleBaseRange(String baseRange) {
		// weird bug with a user getting whitespace in their string somehow.
		baseRange.replaceAll("\\s+","");
		// another weird bug for non-ascii characters.  Reject any character that isn't a digit or "-"
		baseRange=sanitizeString(baseRange);
		String [] split = baseRange.split("-");
		if (split.length!=2)
			throw new IllegalArgumentException("Unable to split input base range into a start and end location: " + baseRange+ " Please format your base range input properly.");
		BaseRange r = new BaseRange(Integer.parseInt(split[0]), Integer.parseInt(split[1]));
		return r;
	}

	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}

	public static String getSequenceForBaseRange (final List<BaseRange> baseRange, final String sequence) {
		StringBuilder result = new StringBuilder();
		for (BaseRange b: baseRange) {
			String s = sequence.substring(b.start-1, b.end);
			result.append(s);
		}
		return result.toString();
	}

	public static byte []  getBytesForBaseRange (final List<BaseRange> baseRange, final byte [] sequence) {
		// this is the expected size of the output sequence.
		int rangeSize = BaseRange.getTotalRangeSize(baseRange);
		byte [] result = new byte [rangeSize];

		int destPos=0;
		for (BaseRange b: baseRange) {
			byte [] subList = Arrays.copyOfRange(sequence, b.start-1, b.end);
			System.arraycopy(subList, 0, result, destPos, subList.length);
			destPos+=subList.length;
		}
		return result;
	}

	public static String sanitizeString (final String input) {
		StringBuilder result = new StringBuilder();
		for(char val : input.toCharArray())
			// 0-9 or "-"
		    if((val >=48 && val<=57) || val==45) result.append(val);
		return (result.toString());
	}

	/**
	 * Create a list of base ranges that are the negative space of the set passed in - they occupy all the bases that aren't in the passed in list.
	 *
	 * @param baseRanges A list of base ranges
	 * @param sequenceLength The total length of the sequence.  A BaseRange object will be generated that covers the gap
	 * between the end of the last BaseRange handed in and this position, if the sequenceLength > last BaseRange.end.
	 * @return
	 */
	public static List<BaseRange> invert (final List<BaseRange> baseRanges, final int sequenceLength) {
		// logic ripped out of IntervalList.invert().
		// wish I could use that code instead, but I don't want to be tied to a fake SAMFileHeader..
		List<BaseRange> result = new ArrayList<BaseRange>();
		Integer lastCoveredPosition = 0; //start at beginning of sequence
        //iterate over list of intervals

            for (final BaseRange i : baseRanges) {
                if (i.getStart() > lastCoveredPosition + 1) //if there's space between the last interval and the current one, add an interval between them
                	result.add(new BaseRange(lastCoveredPosition + 1, i.getStart()-1));
                lastCoveredPosition = i.getEnd(); //update the last covered position
            }
        //finally, if there's room between the last covered position and the end of the sequence, add an interval
        if (sequenceLength > lastCoveredPosition) //if there's space between the last interval and the next
            // one, add an interval. This also covers the case that there are no intervals in the ListMap for a contig.
            result.add(new BaseRange(lastCoveredPosition + 1, sequenceLength));
        return result;
    }

	public static int getTotalRangeSize(final List<BaseRange> ranges) {
		int result =0;
		for (BaseRange br: ranges) {
			int l = (br.end-br.start)+1;
			result+=l;
		}
		return result;
	}

	@Override
	public String toString () {
		return new String (this.start+"-"+this.end);
	}
}





