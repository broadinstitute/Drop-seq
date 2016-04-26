package org.broadinstitute.dropseqrna.utils;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;

import java.util.Iterator;

public class VariantContextSingletonFilter extends FilteredIterator<VariantContext>{

	private final boolean hetVarOnly;

	/**
	 * Filters out sites that have > 1 sample with an alternate allele
	 * @param underlyingIterator
	 * @param hetVarOnly True if only heterozygous sites are allowed for a variant site.  If true, observing a sample with the alternate allele will reject the record.
	 * If false then either the heterozygous or homozygous alternate will count as a variant observation.
	 * Otherwise count heterozygous and hom_var sites.
	 */
	public VariantContextSingletonFilter(final Iterator<VariantContext> underlyingIterator, final boolean hetVarOnly) {
		super(underlyingIterator);
		this.hetVarOnly=hetVarOnly;

	}

	@Override
	/**
	 * For each variant context, look at the genotypes of the samples, and count the number of samples that have
	 * an alternate allele.  If that number of samples!=1, return true to filter this record.
	 * @param rec
	 * @return
	 */
	public boolean filterOut(final VariantContext rec) {
		GenotypesContext gc = rec.getGenotypes();
		Iterator<Genotype> iter = gc.iterator();
		int count=0;
		while (iter.hasNext()) {
			Genotype g = iter.next();
			// boolean isHet = g.isHet();
			// boolean homRef = g.isHomRef();

			if (hetVarOnly & g.isHomVar())
				return true; // filter when het only and observe hom_var.
			if (g.isHet() || g.isHomVar())
				count++;
		}
		if (count==1) return false;
		return true;
	}

}
