/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.dropseqrna.utils.alignmentcomparison;

import java.util.List;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;

/**
 * Holds two iterators and outputs data when both iterators have the same query
 * name. When one iterator is "ahead" of the other, the only the lagging
 * iterator is incremented.
 *
 * @author nemesh
 *
 */
public class QueryNameJointIterator {

	private static final Log log = Log.getInstance(QueryNameJointIterator.class);

	private final PeekableIterator<List<SAMRecord>> iterOne;
	private final PeekableIterator<List<SAMRecord>> iterTwo;
	private final SAMRecordQueryNameComparator comp;

	private JointResult next = null;

	public QueryNameJointIterator(final PeekableIterator<List<SAMRecord>> iterOne, final PeekableIterator<List<SAMRecord>> iterTwo) {
		this.comp = new SAMRecordQueryNameComparator();
		this.iterOne = iterOne;
		this.iterTwo = iterTwo;
		getNextSet();
	}



	public boolean hasNext() {
		if (this.next == null) getNextSet(); // iterates until you have a result, or you're out of results.
		return this.next != null;
	}

	/**
	 * Get the next result.
	 *
	 * @return
	 */
	public JointResult next() {
		if (this.next == null)
			getNextSet(); // iterates until you have a result, or you're out of
							// results.
		// the result you'll return
		JointResult result = this.next;
		// since you are handing out a result, the cached result is null.
		this.next = null;
		return result;
	}

	private void getNextSet() {

		while (iterOne.hasNext() && iterTwo.hasNext()) {

			List<SAMRecord> r1List = iterOne.peek();
			List<SAMRecord> r2List = iterTwo.peek();
			// if (r1List.get(0).getReadName().equals(anObject))
			// only have to compare the first record.
			int cmp = comp.fileOrderCompare(r1List.get(0), r2List.get(0));

			if (cmp < 0)
				r1List = iterOne.next();
			else if (cmp > 0)
				r2List = iterTwo.peek();
			else if (cmp == 0) {
				// do some real work.
				// grab the next record and process it.
				r1List = iterOne.next();
				r2List = iterTwo.next();
				JointResult jr = new JointResult(r1List, r2List);
				this.next=jr;
				break;
			}
		}
	}

	public class JointResult {

		private final List<SAMRecord> one;
		private final List<SAMRecord> two;

		public JointResult(final List<SAMRecord> one, final List<SAMRecord> two) {
			if (one.isEmpty() || two.isEmpty())
				log.warn("One record list empty, this shouldn't happen...");
			if (!one.get(0).getReadName().equals(two.get(0).getReadName()))
				log.error("Read names don't match for lists....");
			this.one = one;
			this.two = two;
		}

		public List<SAMRecord> getOne() {
			return this.one;
		}

		public List<SAMRecord> getTwo() {
			return this.two;
		}

	}
}
