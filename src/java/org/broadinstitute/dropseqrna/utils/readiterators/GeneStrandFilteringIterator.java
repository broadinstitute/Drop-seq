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
package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.Iterator;

import org.broadinstitute.dropseqrna.barnyard.Utils;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;

import htsjdk.samtools.SAMRecord;

/**
 * Filters reads where the read strand and gene strand disagree.
 * @author nemesh
 *
 */
public class GeneStrandFilteringIterator extends FilteredIterator<SAMRecord>{

	private final String strandTag;

	protected GeneStrandFilteringIterator(final Iterator<SAMRecord> underlyingIterator, final String strandTag) {
		super(underlyingIterator);
		this.strandTag=strandTag;
	}

	@Override
	/**
	 * For a given record, if the gene strand tag does not match the read strand tag, reject the record.
	 * If there are multiple gene strand tags, they must all agree with the read strand tag.
	 * @param rec
	 * @return
	 */
	public boolean filterOut(final SAMRecord rec) {
		/*
		if (rec.getReadName().equals("HN7TNBGXX:4:11609:2753:3252"))
			System.out.println("STOP");
		*/
		String geneStrand = rec.getStringAttribute(strandTag);
		if (geneStrand==null) return false;

		String [] strands = geneStrand.split(",");
		String readStrandString = Utils.strandToString(!rec.getReadNegativeStrandFlag());
		for (String s: strands)
			if (!s.equals(readStrandString)) return true;
		return false;
	}


}
