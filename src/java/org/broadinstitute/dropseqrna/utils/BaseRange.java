package org.broadinstitute.dropseqrna.utils;

import java.util.ArrayList;
import java.util.List;

public class BaseRange {
	
	private int start;
	private int end;

	public BaseRange(int start, int end) {
		this.start= start;
		this.end=end;
	}
	
	public int getSize () {
		int r = this.getEnd()-this.getStart()+1;
		return (r);
	}
	
	public static int getTotalRangeSize (String baseRange) {
		List<BaseRange> result = parseBaseRange(baseRange);
		int totalSize=0;
		for (BaseRange b: result) {
			int s = b.getSize();
			totalSize+=s;
		}
		return totalSize;
	}
	
	public static List<BaseRange> parseBaseRange(String baseRange) {
		String [] split = baseRange.split(":");
		List<BaseRange> result = new ArrayList<BaseRange>(split.length);
		for (String s: split) {
			BaseRange r = parseSingleBaseRange(s);
			result.add(r);
		}
		return result;
	}

	public static BaseRange parseSingleBaseRange(String baseRange) {
		String [] split = baseRange.split("-");
		BaseRange r = new BaseRange(Integer.parseInt(split[0]), Integer.parseInt(split[1]));
		return r;
	}

	public int getStart() {
		return start;
	}

	public int getEnd() {
		return end;
	}
	
	public static String getSequenceForBaseRange (List<BaseRange> baseRange, String sequence) {
		StringBuilder result = new StringBuilder();
		for (BaseRange b: baseRange) {
			String s = sequence.substring(b.start-1, b.end);
			result.append(s);
		}
		return result.toString();		
	}
	

}
