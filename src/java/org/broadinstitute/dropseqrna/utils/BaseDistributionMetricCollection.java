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

import htsjdk.samtools.util.IOUtil;
import org.apache.commons.lang.StringUtils;
import org.broadinstitute.dropseqrna.TranscriptomeException;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class BaseDistributionMetricCollection implements Serializable {

	/**
	 *
	 */
	private static final long serialVersionUID = -55813371198085699L;

	private Map<Integer, BaseDistributionMetric> collection = null;

	public BaseDistributionMetricCollection() {
		collection = new HashMap<>();
	}

    public void addBase (final char base, final int position, final int numTimes) {
        BaseDistributionMetric m = this.collection.get(position);
        if (m==null) {
            m = new BaseDistributionMetric();
            this.collection.put(position, m);
        }
        m.addBase(base, numTimes);
    }

	public void addBase (final char base, final int position) {
        addBase(base, position, 1);
	}

	public void addBases (final String bases) {
		char [] b = bases.toCharArray();
		for (int i=0; i<b.length; i++)
			addBase(b[i], i);
	}

	public void addBases (final char [] bases) {
		for (int i=0; i<bases.length; i++)
			addBase(bases[i], i);
	}

	public void addBases (final byte [] bases) {
		for (int i=0; i<bases.length; i++) {
			char b = (char) bases[i];
			addBase(b, i);
		}
	}

    public void mergeMetricCollections(BaseDistributionMetricCollection otherCollection) {
        if (otherCollection == null)
            return;

        for (int position : otherCollection.getPositions()) {
            BaseDistributionMetric metric = getDistributionAtPosition(position);
            if (metric == null) {
                this.collection.put(position, otherCollection.getDistributionAtPosition(position));
            } else {
                metric.mergeMetrics(otherCollection.getDistributionAtPosition(position));
            }
        }
    }

	public BaseDistributionMetric getDistributionAtPosition (final int position) {
		return this.collection.get(position);
	}

	public List<Integer> getPositions () {
		List<Integer> r = new ArrayList<>(this.collection.keySet());
		Collections.sort(r);
		return (r);
	}

    public void writeOutput (final File output) {
        BufferedWriter writer = OutputWriterUtil.getWriter(output);
        String [] header = {"position", "A", "C","G", "T", "N"};
        String h = StringUtils.join(header, "\t");
        OutputWriterUtil.writeResult(h, writer);

        List<Integer> sortedKeys = getPositions();
        for (Integer key : sortedKeys) {
            BaseDistributionMetric brd = getDistributionAtPosition(key);

            String [] l = {key + "",
                    brd.getCount(Bases.A.getBase())+"", brd.getCount(Bases.C.getBase())+"",
                    brd.getCount(Bases.G.getBase())+"",brd.getCount(Bases.T.getBase())+"",
                    brd.getCount(Bases.N.getBase())+""};
            String line = StringUtils.join(l, "\t");
            OutputWriterUtil.writeResult(line, writer);
        }
        OutputWriterUtil.closeWriter(writer);
    }

    public static BaseDistributionMetricCollection readBaseDistribution(final File reportFile) {
        BaseDistributionMetricCollection metricCollection = new BaseDistributionMetricCollection();

        String[] columns = null;
        char[] bases = null;
        try {
            BufferedReader input = IOUtil.openFileForBufferedReading(reportFile);
            try  {
                String line;
                while ((line = input.readLine()) != null) {
                    line=line.trim();
                    if (columns == null && line.startsWith("position")) {
                        columns = line.split("\t");
                        if (columns.length != 6) {
                            throw new RuntimeException("Report header for " + reportFile.getAbsolutePath() + " contains " + columns.length + " columns");
                        }
                        bases = new char[columns.length];
                        for (int idx=1; idx<columns.length; idx++) {
                            bases[idx] = columns[idx].charAt(0);
                        }
                        continue;
                    } else if (columns == null) {
                        throw new RuntimeException("The first line in " + reportFile.getAbsolutePath() + " does not look like a header");
                    }
                    String[] strLine = line.split("\t");
                    int position = Integer.parseInt(strLine[0]);
                    for (int idx=1; idx<6; idx++) {
                        int count = Integer.parseInt(strLine[idx]);
                        metricCollection.addBase(bases[idx], position, count);
                    }
                }
            } finally {
                input.close();
            }
        } catch (IOException ex) {
            throw new TranscriptomeException("Error reading the file: " + reportFile.getAbsolutePath());
        }

        return (metricCollection);
    }

    @Override
	public String toString () {
		return this.collection.toString();
	}

    public boolean distributionsEqual(final BaseDistributionMetricCollection other) {
        if (other.getPositions().size() != getPositions().size())
            return false;

        for (int position : getPositions()) {
            if (!getDistributionAtPosition(position).distributionsEqual(other.getDistributionAtPosition(position)))
                return false;
        }

        return true;
    }

}
