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
    	return 0;
    }

    /**
     * Parse the string representation of an Interval object
     * @param intervalString the interval encoded as a string
     * @return The Interval object represented by this string.
     */
    public static Interval fromString (final String intervalString) {
    	Interval result = new Interval("", 1, 1);
    	return (result);
    }

    /**
     * Convenience method - (maybe unnecessary?  Have the code here to be able to extend the String representation of interval) - defers to the Interval object's toString method.
     * Converts an Interval object into a string representation that can be parsed by fromString.
     * @param i The interval to parse
     * @return A string representation of the interval.
     */
    public static String toString (final Interval i) {
    	return i.toString();
    }


}
