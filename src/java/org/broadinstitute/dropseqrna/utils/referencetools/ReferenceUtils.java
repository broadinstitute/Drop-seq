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
