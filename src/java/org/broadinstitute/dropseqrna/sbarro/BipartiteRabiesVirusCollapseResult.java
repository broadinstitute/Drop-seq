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

package org.broadinstitute.dropseqrna.sbarro;

/**
 * For a single rabies sequence ("child"), holds a potential larger "parent" sequence, the edit distance, and the UMI sharing result.
 *
 * @author nemesh
 *
 */
public class BipartiteRabiesVirusCollapseResult {

	private final String childBarcode;
	private final String parentBarcode;
	private final int editDistanceLeft;
	private final int editDistanceRight;
	private Double umiSharing;
	private Integer umisChild;
	private Integer umisParent;

	public BipartiteRabiesVirusCollapseResult(final String childBarcode, final String parentBarcode,
			final int editDistanceLeft, final int editDistanceRight) {
		this.childBarcode=childBarcode;
		this.parentBarcode=parentBarcode;
		this.editDistanceLeft=editDistanceLeft;
		this.editDistanceRight=editDistanceRight;
	}

	public BipartiteRabiesVirusCollapseResult(final String childBarcode, final String parentBarcode,
			final int editDistanceLeft, final int editDistanceRight, final double umiSharing) {
		this(childBarcode, parentBarcode, editDistanceLeft, editDistanceRight);
		this.umiSharing=umiSharing;
	}

	public String getChildBarcode() {
		return childBarcode;
	}

	public String getParentBarcode() {
		return parentBarcode;
	}

	public Double getUmiSharing() {
		return umiSharing;
	}

	public void setUmiSharing(final Double value) {
		this.umiSharing=value;
	}

	public int getEditDistanceLeft() {
		return editDistanceLeft;
	}

	public int getEditDistanceRight() {
		return editDistanceRight;
	}

	public Integer getUmisChild() {
		return umisChild;
	}

	public void setUmisChild(final Integer value) {
		this.umisChild=value;
	}

	public Integer getUmisParent() {
		return umisParent;
	}

	public void setUmisParent(final Integer value) {
		this.umisParent=value;
	}

	@Override
	public String toString () {
		StringBuilder b = new StringBuilder();
		b.append("For child barcode [" +this.getChildBarcode() +"]  Best parent [" + this.getParentBarcode()+"]");
		b.append(" left ED [" + this.getEditDistanceLeft() +"] right ED [" + this.getEditDistanceRight()+"]");
		if (this.getUmiSharing()!=null) b.append(" UMI Sharing [" +  Double.toString(this.getUmiSharing())+"]");
		if (this.getUmisChild()!=null) b.append(" UMIs child [" + this.getUmisChild()+"]");
		if (this.getUmisParent()!=null) b.append(" UMIs parent [" + this.getUmisParent()+"]");
		return b.toString();
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((childBarcode == null) ? 0 : childBarcode.hashCode());
		result = prime * result + editDistanceLeft;
		result = prime * result + editDistanceRight;
		result = prime * result + ((parentBarcode == null) ? 0 : parentBarcode.hashCode());
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
		BipartiteRabiesVirusCollapseResult other = (BipartiteRabiesVirusCollapseResult) obj;
		if (childBarcode == null) {
			if (other.childBarcode != null)
				return false;
		} else if (!childBarcode.equals(other.childBarcode))
			return false;
		if (editDistanceLeft != other.editDistanceLeft)
			return false;
		if (editDistanceRight != other.editDistanceRight)
			return false;
		if (parentBarcode == null) {
			if (other.parentBarcode != null)
				return false;
		} else if (!parentBarcode.equals(other.parentBarcode))
			return false;
		return true;
	}


}
