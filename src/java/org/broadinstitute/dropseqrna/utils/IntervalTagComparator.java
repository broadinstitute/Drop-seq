package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMTagUtil;
import htsjdk.samtools.util.Interval;

import java.util.Comparator;

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

    @Override
    public int compare(final SAMRecord rec1, final SAMRecord rec2) {
    	Interval i1 = getIntervalForTag(rec1);
    	Interval i2 = getIntervalForTag(rec2);
    	return (compare(i1,i2, this.dict));
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
    	// at most there are 5 fields.  This takes care of having any characters in names.
    	String [] s = intervalString.split(":", 4);
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
    	start=Integer.parseInt(posArray[0]);
    	// set the end = the start, this changes if there's a second position in the field.
    	end=Integer.parseInt(posArray[0]);
    	if (posArray.length==2)
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
     * Instead, use ":" as a delimiter for all fields except the range fields start-end.
     * @param i The interval to parse
     * @return A string representation of the interval.
     */
    public static String toString (final Interval i) {
    	// chr1:1-10	+	foo
    	String result = i.getContig() + ":" + i.getStart() + "-" + i.getEnd() + ":" + (i.isNegativeStrand() ? '-' : '+') + ":" + ((null == i.getName()) ? '.' : i.getName());
    	return result;
    }


}
