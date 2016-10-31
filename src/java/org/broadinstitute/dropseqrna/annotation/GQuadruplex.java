package org.broadinstitute.dropseqrna.annotation;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class GQuadruplex {

	private Interval interval;
	private String [] sequence;
	private static String patternString="([gG]{3,})([acgtACGT]{1,7})([gG]{3,})([acgtACGT]{1,7})([gG]{3,})([acgtACGT]{1,7})([gG]{3,})";
	private static Pattern pattern=Pattern.compile(patternString);

	public GQuadruplex (final Interval i, final String G1, final String L1, final String G2, final String L2, final String G3, final String L3, final String G4) {
		this.interval=i;
		sequence=new String [7];
		sequence[0]=G1;
		sequence[1]=L1;
		sequence[2]=G2;
		sequence[3]=L2;
		sequence[4]=G3;
		sequence[5]=L3;
		sequence[6]=G4;
	}

	public Interval getMatchInterval() {
		return interval;
	}


	public String getG1 () {
		return sequence[0];
	}

	public String getG2 () {
		return sequence[2];
	}

	public String getG3 () {
		return sequence[4];
	}

	public String getG4 () {
		return sequence[6];
	}

	public String getL1 () {
		return sequence[1];
	}

	public String getL2 () {
		return sequence[3];
	}

	public String getL3 () {
		return sequence[5];
	}

	public String getSequence() {
		return StringUtil.join("", this.sequence);
	}


	public static List<GQuadruplex> find (final String seqName, final String seq) {
		Matcher matcher = pattern.matcher(seq);
		List<GQuadruplex> result = new ArrayList<GQuadruplex>();
		while (matcher.find()) {
			int start = matcher.start();
			int end = matcher.end();
			// Make the start and end 1 based.
			Interval i = new Interval(seqName, start+1, end);
			GQuadruplex r = new GQuadruplex(i, matcher.group(1), matcher.group(2), matcher.group(3), matcher.group(4), matcher.group(5), matcher.group(6), matcher.group(7));
			result.add(r);
		}
		return (result);
	}

	@Override
	public String toString () {
		StringBuilder b = new StringBuilder();
		b.append(this.interval.getContig()+":"+this.interval.getStart()+"-" +this.interval.getEnd());
		b.append(" "+ this.getSequence()+" ");
		b.append(" G1 [" + getG1()+"]");
		b.append(" L1 [" + getL1() +"]");

		b.append(" G2 [" + getG2()+"]");
		b.append(" L2 [" + getL2() +"]");

		b.append(" G3 [" + getG3()+"]");
		b.append(" L3 [" + getL3() +"]");

		b.append(" G4 [" + getG4()+"]");
		return b.toString();
	}
}
