/*
 * MIT License
 *
 * Copyright 2019 Broad Institute
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
package org.broadinstitute.dropseqrna.metrics;

import htsjdk.samtools.metrics.MetricBase;

public class UmiSharingMetrics
        extends MetricBase {
    /** collapse tag value */
    public String PARENT;

    /** uncollapsed tag value */
    public String CHILD;

    /** number of distinct UMIs (or other COUNT_TAG value) in parent reads */
    public int NUM_PARENT;

    /** number of distinct UMIs (or other COUNT_TAG value) in this group of child reads */
    public int NUM_CHILD;

    /** number of distinct child UMIs that match a parent UMI, possibly with edit distance */
    public int NUM_SHARED;

    /** NUM_SHARED/NUM_CHILD */
    public double FRAC_SHARED;

    public String toString () {
    	return "Parent [" + this.PARENT +"] Child ["+this.CHILD+"] NUM_PARENT ["+this.NUM_PARENT+"] NUM_CHILD [" + this.NUM_CHILD +"] NUM_SHARED ["+this.NUM_SHARED+"] FRAC_SHARED [" + this.FRAC_SHARED+"]";
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = super.hashCode();
		result = prime * result + ((CHILD == null) ? 0 : CHILD.hashCode());
		result = prime * result + NUM_CHILD;
		result = prime * result + NUM_PARENT;
		result = prime * result + NUM_SHARED;
		result = prime * result + ((PARENT == null) ? 0 : PARENT.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (!super.equals(obj))
			return false;
		if (getClass() != obj.getClass())
			return false;
		UmiSharingMetrics other = (UmiSharingMetrics) obj;
		if (CHILD == null) {
			if (other.CHILD != null)
				return false;
		} else if (!CHILD.equals(other.CHILD))
			return false;
		if (NUM_CHILD != other.NUM_CHILD)
			return false;
		if (NUM_PARENT != other.NUM_PARENT)
			return false;
		if (NUM_SHARED != other.NUM_SHARED)
			return false;
		if (PARENT == null) {
			if (other.PARENT != null)
				return false;
		} else if (!PARENT.equals(other.PARENT))
			return false;
		return true;
	}

	
    
    
}
