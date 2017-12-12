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

import java.util.Arrays;
import java.util.List;

import com.google.common.collect.ComparisonChain;
import com.google.common.collect.Ordering;

public class ContigResult implements Comparable<ContigResult> {

	/**
	 * Only populate on of newContig or newContigs.
	 * Don't incur the extra memory load of storing a list of contigs containing a single element.
	 */
	private final String oldContig;
	private String newContig;
	private List<String> newContigs;
	private final boolean newReadMapsUniquely;

	public ContigResult (final String oldContig, final String newContig, final boolean newReadMapsUniquely) {
		this.oldContig=oldContig;
		this.newContig=newContig;
		this.newContigs=null;
		this.newReadMapsUniquely=newReadMapsUniquely;
	}

	public ContigResult (final String oldContig, final List<String> newContigs, final boolean newReadMapsUniquely) {
		this.oldContig=oldContig;
		this.newContig=null;
		this.newContigs=newContigs;
		this.newReadMapsUniquely=newReadMapsUniquely;
	}

	public String getOldContig() {
		return oldContig;
	}

	public List<String> getNewContigs() {
		if (this.newContigs!=null) return this.newContigs;
		return Arrays.asList(this.newContig);
	}

	public boolean isNewReadMapsUniquely() {
		return newReadMapsUniquely;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((newContig == null) ? 0 : newContig.hashCode());
		result = prime * result + ((newContigs == null) ? 0 : newContigs.hashCode());
		result = prime * result + (newReadMapsUniquely ? 1231 : 1237);
		result = prime * result + ((oldContig == null) ? 0 : oldContig.hashCode());
		return result;
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		ContigResult other = (ContigResult) obj;
		if (newContig == null) {
			if (other.newContig != null)
				return false;
		} else if (!newContig.equals(other.newContig))
			return false;
		if (newContigs == null) {
			if (other.newContigs != null)
				return false;
		} else if (!newContigs.equals(other.newContigs))
			return false;
		if (newReadMapsUniquely != other.newReadMapsUniquely)
			return false;
		if (oldContig == null) {
			if (other.oldContig != null)
				return false;
		} else if (!oldContig.equals(other.oldContig))
			return false;
		return true;
	}

	@Override
	public int compareTo(final ContigResult o) {
		ComparisonChain.start().compare(this.oldContig, o.getOldContig())
		.compare(this.getNewContigs(), o.getNewContigs(), Ordering.<String>natural().lexicographical())
		.compareTrueFirst(this.newReadMapsUniquely, o.newReadMapsUniquely)
		.result();
		return 0;
	}

	@Override
	public String toString () {
		return "original mapping[" + this.oldContig + "] new mapping(s) [" + this.getNewContigs().toString() +" ] new read unique mapping [" + this.newReadMapsUniquely +"]";
	}



}
