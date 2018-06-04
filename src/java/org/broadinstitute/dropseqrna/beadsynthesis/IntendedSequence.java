/*
 * MIT License
 *
 * Copyright 2018 Broad Institute
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

package org.broadinstitute.dropseqrna.beadsynthesis;

import org.apache.commons.lang.StringUtils;

import java.util.Collections;
import java.util.List;

public class IntendedSequence {

	private final String intendedSequence;
	private final List<String> relatedSequences;
	private final Character deletedBase;
	private final Integer deletedBasePos;
	private final Integer intendedSequenceUMIs;
	private final Double medianRelatedSequenceUMIs;
	private final Integer totalRelatedUMIs;
	private final Double intendedSequenceUMIBias;
	private final Double relatedMedianUMIBias;

	// Intended sequence, related sequences, deleted base, deleted base position, rate, intended base, intended UMis, median neighbor UMIs,
	public IntendedSequence (final String intendedSequence, final List<String> relatedSequences, final Character deletedBase, final Integer deletedBasePos,
			final Integer intendedSequenceUMIs, final Double medianRelatedSequenceUMIs, final Integer totalRelatedUMIs, final Double intendedSequenceUMIBias, final double relatedMedianUMIBias) {
		this.intendedSequence=intendedSequence;
		Collections.sort(relatedSequences);
		this.relatedSequences=relatedSequences;
		this.deletedBase=deletedBase;
		this.deletedBasePos=deletedBasePos;
		this.intendedSequenceUMIs=intendedSequenceUMIs;
		this.medianRelatedSequenceUMIs=medianRelatedSequenceUMIs;
		this.totalRelatedUMIs=totalRelatedUMIs;
		this.intendedSequenceUMIBias=intendedSequenceUMIBias;
		this.relatedMedianUMIBias=relatedMedianUMIBias;
	}

	public String getIntendedSequence() {
		return intendedSequence;
	}

	public List<String> getRelatedSequences() {
		return relatedSequences;
	}

	public Character getDeletedBase() {
		return deletedBase;
	}

	public Integer getDeletedBasePos() {
		return deletedBasePos;
	}

	public Double getDeletionRate() {
		if (this.deletedBase==null) return null;
		double dr=(double) this.totalRelatedUMIs/  (double) (this.intendedSequenceUMIs+this.totalRelatedUMIs);
		return dr;
	}

	public Integer getIntendedSequenceUMIs() {
		return intendedSequenceUMIs;
	}

	public Double getMedianRelatedSequenceUMIs() {
		return medianRelatedSequenceUMIs;
	}

	public Double getIntendedSequenceUMIBias() {
		return intendedSequenceUMIBias;
	}

	public Double getRelatedMedianUMIBias() {
		return relatedMedianUMIBias;
	}

	@Override
	public String toString () {
		StringBuilder b = new StringBuilder();
		b.append("Intended Sequence [" + getNullableString(this.intendedSequence) +"] ");
		b.append("Neighbor Sequences " + StringUtils.join(this.relatedSequences, ":")+ " ");
		if (this.deletedBase!=null) b.append("deleted base [" + this.deletedBase +"] ");
		if (this.deletedBasePos!=null) b.append("deleted base pos [" + this.deletedBasePos +"] ");
		if (this.intendedSequenceUMIs!=null) b.append("intended sequence UMIs [" + this.intendedSequenceUMIs+"] ");
		b.append("Median related UMIs [" + this.medianRelatedSequenceUMIs +"] ");
		b.append("Total related UMIs [" + this.totalRelatedUMIs+"] ");
		if (this.intendedSequenceUMIBias!=null) b.append("Intended UMI Bias ["+ intendedSequenceUMIBias+"] ");
		b.append("Related Median UMI T Bias [" + this.relatedMedianUMIBias+"]");
		return b.toString();
	}

	private String getNullableString (final Object o) {
		if (o==null) return "";
		return o.toString();
	}


}
