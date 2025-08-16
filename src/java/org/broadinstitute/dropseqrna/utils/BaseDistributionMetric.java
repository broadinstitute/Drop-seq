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

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.metrics.MetricBase;
import htsjdk.samtools.util.Log;

public class BaseDistributionMetric extends MetricBase implements Serializable {

	/**
	 *
	 */
	private static final long serialVersionUID = -644875076890522563L;

    private static final Log LOG = Log.getInstance(BaseDistributionMetric.class);

	private Map<Character, Integer> map = null;

	public BaseDistributionMetric() {
		map = new HashMap<>();
		map.put(Bases.A.getBase(), 0);
		map.put(Bases.C.getBase(), 0);
		map.put(Bases.G.getBase(), 0);
		map.put(Bases.T.getBase(), 0);
		map.put(Bases.N.getBase(), 0);
	}

	public BaseDistributionMetric(final int countA, final int countC, final int countG, final int countT, final int countN) {
		map = new HashMap<>();
		map.put(Bases.A.getBase(), countA);
		map.put(Bases.C.getBase(), countC);
		map.put(Bases.G.getBase(), countG);
		map.put(Bases.T.getBase(), countT);
		map.put(Bases.N.getBase(), countN);
	}

    void addBase(final Character base, final int numTimes, final ValidationStringency validationStringency) {
        Integer count = map.get(base);
        if (count != null) {
            count += numTimes;
            map.put(base, count);
        } else {
            switch (validationStringency) {
                case STRICT:
                    throw new IllegalArgumentException("Base '" + base + "' not found in BaseDistributionMetric map.");
                case LENIENT:
                    LOG.warn("Base '" + base + "' not found in BaseDistributionMetric map.");
                    break;
                case SILENT:
                    break;
            }
            addBase(Bases.N.getBase(), numTimes, validationStringency);
        }
    }

    void addBase(final Character base, final int numTimes) {
        addBase(base, numTimes, ValidationStringency.DEFAULT_STRINGENCY);
    }

    void addBase(final Character base) {
        addBase(base, 1);
	}

    public int getCount(final Character base) {
		return (map.get(base));
	}

    public int getTotalCount () {
		int count=0;
		for (Bases b: Bases.values()) {
			char bb = b.getBase();
			count+=map.get(bb);
		}
		return count;
	}

    public void mergeMetrics(BaseDistributionMetric otherMetric) {
        for (Bases base : Bases.values()) {
            addBase(base.getBase(), otherMetric.getCount(base.getBase()));
        }
    }

    @Override
	public String toString() {
		return this.map.toString();
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = super.hashCode();
		result = prime * result + ((map == null) ? 0 : map.hashCode());
		return result;
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (!super.equals(obj))
			return false;
		if (getClass() != obj.getClass())
			return false;
		BaseDistributionMetric other = (BaseDistributionMetric) obj;
		if (map == null) {
			if (other.map != null)
				return false;
		} else if (!map.equals(other.map))
			return false;
		return true;
	}

    public boolean distributionsEqual(final BaseDistributionMetric otherMetric) {
        for (Bases base : Bases.values()) {
            if (otherMetric.getCount(base.getBase()) != getCount(base.getBase()))
                return false;
        }

        return true;
    }
}
