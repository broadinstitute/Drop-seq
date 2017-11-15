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

import java.util.Comparator;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTagUtil;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;

/**
 * When assigning intervals as tags to a BAM, the interval information is stored as a string.
 * This Comparator allows sorting of various flavors of interval String encodings,
 * by converting the String version of the interval back to a proper Interval object,
 * then comparing intervals.
 * If a sequence dictionary is supplied, then the use the sequenceIndex integer to sort the objects
 * instead of the String representation of the sequence.
 * @author nemesh
 *
 */
public class IntervalTagComparator implements Comparator<SAMRecord> {


	private final short tag;
	private final SAMSequenceDictionary dict;
	private static final Pattern colonPattern = Pattern.compile(":");
	private static final Pattern dashPattern = Pattern.compile("-");

	public static final char ENCODE_DELIMITER='|';


	/**
	 * Sorts the string representations of intervals
	 * @param tag
	 */
    public IntervalTagComparator(final String tag) {
        this.tag = SAMTagUtil.getSingleton().makeBinaryTag(tag);
        this.dict=null;
    }

    public IntervalTagComparator (final String tag, final SAMSequenceDictionary dict) {
    	this.tag = SAMTagUtil.getSingleton().makeBinaryTag(tag);
    	this.dict=dict;
    }


    public int compareOld(final SAMRecord rec1, final SAMRecord rec2) {
    	Interval i1 = getIntervalForTag(rec1);
    	Interval i2 = getIntervalForTag(rec2);
    	return (compare(i1,i2, this.dict));
    }

    @Override
    public int compare(final SAMRecord rec1, final SAMRecord rec2) {
    	// instead of parsing the interval fully, parse chromosome, then start/end.
    	String [] rec1Split = getFirstSplit(rec1);
    	String [] rec2Split = getFirstSplit(rec2);

    	String rec1Contig=rec1Split[0];
    	String rec2Contig=rec2Split[0];

    	int result = 0;
    	// compare the contig.
    	if (dict!=null) {
    		int seqIdx1 = dict.getSequenceIndex(rec1Contig);
    		int seqIdx2 = dict.getSequenceIndex(rec2Contig);
    		result = seqIdx1 - seqIdx2;
    	} else
			result = rec1Contig.compareTo(rec2Contig);

    	// if the contig isn't enough, parse out the start/end.
    	if (result == 0) {
    		String pos1= rec1Split[1];
    		String pos2= rec2Split[1];

    		int [] posArray1 = parsePosition(pos1);
    		int [] posArray2 = parsePosition(pos2);


            if (posArray1[0] == posArray2[0])
				result = posArray1[1] - posArray2[1];
			else
				result = posArray1[0] - posArray2[0];

            // added bonus, sort on interval names to tie break, if both intervals have names.
            if (result==0) {
            	String n1 = rec1Split[3];
            	String n2 = rec2Split[3];
            	if (n1!=null && n2!=null)
					result = n1.compareTo(n2);
            }
        }

    	return result;

    }

    private String [] getFirstSplit(final SAMRecord rec) {

    	Object strIntervalRec = rec.getAttribute(tag);
    	if (!(strIntervalRec instanceof String))
			throw new IllegalArgumentException(SAMTagUtil.getSingleton().makeStringTag(this.tag) + " does not have a String value");
    	String intervalString = (String) strIntervalRec;
    	String [] result = new String [4];
    	StringUtil.splitConcatenateExcessTokens(intervalString, result, ENCODE_DELIMITER);
    	return (result);
    }

    // parses the start:end position
    private int [] parsePosition (final String posString) {
    	String [] posArray = new String [2];
    	StringUtil.split(posString, posArray, '-');
    	int start=Integer.parseInt(posArray[0]);
    	// set the end = the start, this changes if there's a second position in the field.
    	int end=Integer.parseInt(posArray[0]);
    	if (posArray[1]!=null)
			end = Integer.parseInt(posArray[1]);

    	int [] result = {start,end};
    	return result;

    }

    public static int compare (final Interval i1, final Interval i2, final SAMSequenceDictionary dict) {
    	int result = 0;
    	// if there's a sequence dictionary, compare on the index of the sequence instead of the contig name.
    	if (dict!=null) {
    		int seqIdx1 = dict.getSequenceIndex(i1.getContig());
    		int seqIdx2 = dict.getSequenceIndex(i2.getContig());
    		result = seqIdx1 - seqIdx2;
    	} else
			result = i1.getContig().compareTo(i2.getContig());
    	// same as Interval.compareTo
    	if (result == 0) {
            if (i1.getStart() == i2.getStart())
				result = i1.getEnd() - i2.getEnd();
			else
				result = i1.getStart() - i2.getStart();
            // added bonus, sort on interval names to tie break, if both intervals have names.
            if (result==0) {
            	String n1 = i1.getName();
            	String n2 = i2.getName();
            	if (n1!=null && n2!=null)
					result = n1.compareTo(n2);
            }
        }
    	return result;
    }

    /**
     * Converts the tag to an interval object if that's possible.
     * @param rec
     * @return
     */
    private Interval getIntervalForTag (final SAMRecord rec) {
    	Object strIntervalRec = rec.getAttribute(tag);
    	if (!(strIntervalRec instanceof String))
			throw new IllegalArgumentException(SAMTagUtil.getSingleton().makeStringTag(this.tag) + " does not have a String value");
    	Interval result = fromString((String) strIntervalRec );
    	return (result);
    }

    /**
     * Parse the string representation of an Interval object
     * @param intervalString the interval encoded as a string
     * @return The Interval object represented by this string.
     */
    public static Interval fromString (final String intervalString) {

    	String [] s = new String [4];

    	StringUtil.splitConcatenateExcessTokens(intervalString, s, ENCODE_DELIMITER);

    	if (s.length==1)
			throw new IllegalArgumentException("This is not an interval " + intervalString + ", as there's no '" + ENCODE_DELIMITER+ "' delimiters");

    	String contig=s[0];
    	Integer start=null;
    	Integer end=null;
    	Boolean negativeStrand=null;
    	String name=null;
    	// s[1] can be 1 or 2 in length depending on if end exists.
    	String pos=s[1];

    	// String [] posArray = pos.split("-");

    	// String [] posArray = dashPattern.split(pos);

    	String [] posArray = new String [2];
    	StringUtil.split(pos, posArray, '-');


    	start=Integer.parseInt(posArray[0]);
    	// set the end = the start, this changes if there's a second position in the field.
    	end=Integer.parseInt(posArray[0]);
    	if (posArray[1]!=null)
			end = Integer.parseInt(posArray[1]);
    	if (s[2]!=null){
    		if (s[2].equals("+"))
				negativeStrand=false;
    		if (s[2].equals("-"))
				negativeStrand=true;
    	}
    	if (s[3]!=null)
			name=s[3];

    	// create the interval
    	Interval result = null;
    	if (negativeStrand==null || name==null)
			result = new Interval(contig, start, end);
    	// name null by default, takes care of if name was set or not.
    	if (negativeStrand!=null && name!=null)
			result = new Interval (contig, start, end, negativeStrand, name);

    	return (result);
    }

    /**
     * Parse the string representation of an Interval object
     * @param intervalString the interval encoded as a string
     * @return The Interval object represented by this string.
     */
    public static Interval fromStringOld (final String intervalString) {
    	// at most there are 5 fields.  This takes care of having any characters in names.

    	String [] s = intervalString.split(":", 4);

    	// String [] s = colonPattern.split(intervalString, 4);

    	if (s.length==1)
			throw new IllegalArgumentException("This is not an interval " + intervalString + ", as there's no ':' delimiters");

    	String contig=s[0];
    	Integer start=null;
    	Integer end=null;
    	Boolean negativeStrand=null;
    	String name=null;
    	// s[1] can be 1 or 2 in length depending on if end exists.
    	String pos=s[1];

    	String [] posArray = pos.split("-");

    	// String [] posArray = dashPattern.split(pos);

    	start=Integer.parseInt(posArray[0]);
    	// set the end = the start, this changes if there's a second position in the field.
    	end=Integer.parseInt(posArray[0]);
    	if (posArray[1]!=null)
			end = Integer.parseInt(posArray[1]);
    	if (s.length>2){
    		if (s[2].equals("+"))
				negativeStrand=false;
    		if (s[2].equals("-"))
				negativeStrand=true;
    	}
    	if (s.length>3)
			name=s[3];

    	// create the interval
    	Interval result = null;
    	if (s.length==2)
			result = new Interval(contig, start, end);
    	// name null by default, takes care of if name was set or not.
    	if (s.length>2)
			result = new Interval (contig, start, end, negativeStrand, name);

    	return (result);
    }



    /**
     * Converts an Interval object into a string representation that can be parsed by fromString.
     * Not using the Interval.toString method, because it has tabs, which might really screw up something else...
     * Instead, use ENCODE_DELIMITER as a delimiter for all fields except the range fields start-end.
     * @param i The interval to parse
     * @return A string representation of the interval.
     */
    public static String toString (final Interval i) {
    	// chr1:1-10	+	foo
    	String result = i.getContig() + ENCODE_DELIMITER + i.getStart() + "-" + i.getEnd() + ENCODE_DELIMITER + (i.isNegativeStrand() ? '-' : '+') + ENCODE_DELIMITER + ((null == i.getName()) ? '.' : i.getName());
    	return result;
    }


}
