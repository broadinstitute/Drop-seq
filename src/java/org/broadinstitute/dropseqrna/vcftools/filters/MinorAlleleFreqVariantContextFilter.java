package org.broadinstitute.dropseqrna.vcftools.filters;

import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;

public class MinorAlleleFreqVariantContextFilter extends FilteredIterator <VariantContext> {

	private final Log log = Log.getInstance(MinorAlleleFreqVariantContextFilter.class);

	private final double threshold;
	private final boolean symmetric;
	private final Collection<String> samples;
	private final int GQThreshold;
	private final boolean verbose;

	public MinorAlleleFreqVariantContextFilter(final Iterator<VariantContext> underlyingIterator, final double threshold, final boolean symmetric, final Collection<String> samples, final int GQThreshold) {
		this(underlyingIterator, threshold, symmetric, samples, GQThreshold, false);
	}

	public MinorAlleleFreqVariantContextFilter(final Iterator<VariantContext> underlyingIterator, final double threshold, final boolean symmetric, final Collection<String> samples, final int GQThreshold, final boolean verbose) {
		super(underlyingIterator);
		this.threshold=threshold;
		this.symmetric=symmetric;
		this.samples=samples;
		this.GQThreshold=GQThreshold;
		this.verbose=verbose;
	}

	public MinorAlleleFreqVariantContextFilter(final Iterator<VariantContext> underlyingIterator, final double threshold, final boolean symmetric) {
		this(underlyingIterator, threshold, symmetric, Collections.emptySet(), 30);
	}

	public MinorAlleleFreqVariantContextFilter(final Iterator<VariantContext> underlyingIterator, final double threshold) {
		this(underlyingIterator, threshold, false, Collections.emptySet(),30);
	}

	@Override
	public boolean filterOut(final VariantContext site) {
		double siteMAF = calculateMinorAlleleFrequency(site, this.GQThreshold, this.samples);
		if (siteMAF<threshold) {
			if (verbose) log.info("Rejecting variant allele freq  [" + siteMAF +"] less than threshold [" + threshold +"] " +site.toStringWithoutGenotypes());
			return true;
		}
		if (symmetric && ((1-siteMAF)<threshold)) {
			if (verbose) log.info("Rejecting variant allele freq  [" + siteMAF +"] greater than 1-threshold [" + (1-threshold) +"] " +site.toStringWithoutGenotypes());
			return true;
		}
		return false;
	}

	/**
	 * Calculate minor allele freq from a subset of donors in VCF.
	 * @param site
	 * @param samples.  If empty, use all samples in VCF.
	 * @return
	 */
	public static double calculateMinorAlleleFrequency(final VariantContext site, final int gqThreshold, final Collection <String> samples) {
		if (samples.isEmpty()) return calculateMinorAlleleFrequency(site, gqThreshold);
		ObjectCounter<Integer> c = new ObjectCounter<>();
		for (String s: samples) {
			Genotype g = site.getGenotype(s);
			if (g!=null & g.isCalled() & g.getGQ()>= gqThreshold) {
				if (g.isHomRef()) c.increment(0);
				if (g.isHomVar()) c.increment(2);
				if (g.isHet()) c.increment(1);
			}
		}
		int total = (c.getCountForKey(0)*2)+(c.getCountForKey(1)*2)+(c.getCountForKey(2)*2);
		int altAlleleCounts = c.getCountForKey(1)+(c.getCountForKey(2)*2);
		double maf = (double) altAlleleCounts / (double) total;
		return maf;
	}

	/**
	 * Calculate minor allele freq from all donors in VCF.
	 * @param site
	 * @return
	 */
	public static double calculateMinorAlleleFrequency(final VariantContext site, final int gqThreshold) {
		return calculateMinorAlleleFrequency(site, gqThreshold, site.getSampleNames());
	}

}
