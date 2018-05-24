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
package org.broadinstitute.dropseqrna.cluster;

import java.util.*;
import java.util.regex.Pattern;

public class GeneEnumerator {
    private final Pattern[] regexps;

    private final ArrayList<String> genes = new ArrayList<>();
    private final Map<String, Integer> geneMap = new HashMap<>();

    public GeneEnumerator(final List<String> filteredGeneRegexps) {
        regexps = new Pattern[filteredGeneRegexps.size()];
        for (int i = 0; i < regexps.length; ++i) {
            regexps[i] = Pattern.compile(filteredGeneRegexps.get(i));
        }
    }

    public int getGeneIndex(final String gene) {
        for (final Pattern pattern : regexps) {
            if (pattern.matcher(gene).find()) {
                return -1;
            }
        }
        Integer index = geneMap.get(gene);
        if (index == null) {
            index = genes.size();
            genes.add(gene);
            geneMap.put(gene, index);
        }
        return index;
    }

    public List<String> getGenes() {
        return Collections.unmodifiableList(genes);
    }

    public String getGeneName(final int geneIndex) {
        return genes.get(geneIndex);
    }
}
