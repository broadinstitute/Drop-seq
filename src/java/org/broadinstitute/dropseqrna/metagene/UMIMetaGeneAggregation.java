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
package org.broadinstitute.dropseqrna.metagene;

import org.broadinstitute.dropseqrna.metagene.MetaGene;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;


/**
 * Aggregates meta gene results across many UMIs.
 * Instead of a UMI voting for a set of unique genes, as well as unambiguous and ambiguous meta genes that each gene one count,
 * let each UMI vote ones for each gene and build up evidence.
 *
 *
 * @author nemesh
 *
 */
public class UMIMetaGeneAggregation {

	private ObjectCounter<String> uniqueGenes;
	private ObjectCounter<MetaGene> ambiguousGenes;
	private ObjectCounter<MetaGene> unAmbiguousGenes;

	public UMIMetaGeneAggregation() {
		uniqueGenes = new ObjectCounter<>();
		ambiguousGenes = new ObjectCounter<>();
		unAmbiguousGenes = new ObjectCounter<>();
	}

    public void add(final UMIMetaGeneAggregation otherAggregator) {
        uniqueGenes.increment(otherAggregator.getUniqueGenes());
        ambiguousGenes.increment(otherAggregator.getAmbiguousGenes());
        unAmbiguousGenes.increment(otherAggregator.getUnAmbiguousGenes());
    }

	public void add(final UMIMetaGeneCollection c, final boolean useEditDistance) {
		c.getUnambiguousMetaGenes().forEach(mg -> unAmbiguousGenes.increment(mg));
		c.getAmbiguousMetaGenes().forEach(mg -> ambiguousGenes.increment(mg));
		c.getUniquelyMappedGenes().forEach(g -> uniqueGenes.increment(g));
	}

	public ObjectCounter<String> getUniqueGenes() {
		return uniqueGenes;
	}

	public ObjectCounter<MetaGene> getAmbiguousGenes() {
		return ambiguousGenes;
	}

	public ObjectCounter<MetaGene> getUnAmbiguousGenes() {
		return unAmbiguousGenes;
	}

	/**
	 * Based on the experiment wide data (this object), select a set of meta genes
	 * where the ratio of unambiguous UMIs : unique UMIs is greater than the provided threshold
	 * @param threshold the ratio of unambiguous UMIs to unique UMIs must be greater than this.
	 * @return
	 */
	public Set<MetaGene> getMetaGeneList (final double threshold) {
		Set<MetaGene> result = new HashSet<>();
		List<MetaGene> metaGenes = this.unAmbiguousGenes.getKeysOrderedByCount(true);
		for (MetaGene mg: metaGenes) {
			int mgCount = this.unAmbiguousGenes.getCountForKey(mg);
			Set<String> genes = mg.getGeneNames();
        	int totalUniqueCounts = genes.stream().mapToInt(x->this.getUniqueGenes().getCountForKey(x)).sum();
        	if (totalUniqueCounts==0) totalUniqueCounts=1; // to fix ratios where the denominator is 0.
        	double ratio = (double) mgCount / (double) totalUniqueCounts;
        	if (ratio>=threshold) result.add(mg);

		}
		return (result);
	}
	
	/**
	 * Based on the experiment wide data (this object), select a set of meta genes
	 * where the ratio of unambiguous UMIs : unique UMIs is greater than the provided threshold
	 * @param threshold the ratio of unambiguous UMIs to unique UMIs must be greater than this.
	 * @return A copy of this object, filtered to the subset of metagenes
	 */
	public UMIMetaGeneAggregation filterDataByMetaGenes (Set<MetaGene> metaGenes) {
		UMIMetaGeneAggregation result = new UMIMetaGeneAggregation();
		
		for (MetaGene mg: metaGenes) {
			result.ambiguousGenes.incrementByCount(mg, this.ambiguousGenes.getCountForKey(mg));
			result.unAmbiguousGenes.incrementByCount(mg, this.unAmbiguousGenes.getCountForKey(mg));
        	mg.getGeneNames().stream().forEach(x -> result.uniqueGenes.setCount(x, this.uniqueGenes.getCountForKey(x)));        	
		}
		
		return result;
	}

    @Override
    public boolean equals(Object otherObject) {
        if (this == otherObject)
            return true;
        if (!(otherObject instanceof UMIMetaGeneAggregation))
            return false;

        final UMIMetaGeneAggregation otherAggregation = (UMIMetaGeneAggregation) otherObject;
        return
                Objects.equals(this.uniqueGenes, otherAggregation.uniqueGenes) &&
                Objects.equals(this.ambiguousGenes, otherAggregation.ambiguousGenes) &&
                Objects.equals(this.unAmbiguousGenes, otherAggregation.unAmbiguousGenes);
    }

    public int hashCode() {
        return Objects.hash(uniqueGenes, ambiguousGenes, unAmbiguousGenes);
    }

/**
	 * Probably want a mapping from each gene to meta genes that it's involved with.
	 * Or a way to iterate across unambiguous meta genes and get counts of ambiguous meta gene evidence.
	 */



}
