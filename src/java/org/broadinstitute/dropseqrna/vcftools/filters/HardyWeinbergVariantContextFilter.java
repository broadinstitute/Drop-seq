package org.broadinstitute.dropseqrna.vcftools.filters;

import htsjdk.samtools.util.Log;
import htsjdk.tribble.util.popgen.HardyWeinbergCalculation;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;

import java.util.Iterator;

public class HardyWeinbergVariantContextFilter extends FilteredIterator <VariantContext> {

	private final Log log = Log.getInstance(HardyWeinbergVariantContextFilter.class);
	private final boolean verbose;
	private double threshold;

	public HardyWeinbergVariantContextFilter(final Iterator<VariantContext> underlyingIterator, final double threshold) {
		this(underlyingIterator, threshold, false);
	}
	public HardyWeinbergVariantContextFilter(final Iterator<VariantContext> underlyingIterator, final double threshold, final boolean verbose) {
		super(underlyingIterator);
		this.threshold=threshold;
		this.verbose=verbose;
	}

	@Override
	public boolean filterOut(final VariantContext site) {
		Double pval = getHardyWeinbergPvalue(site);
		if (pval==null) return false; // only reject sites that have some data so you can reject it.
		if (pval<threshold) {
			if (verbose) log.info("Rejecting variant HWE pval [" + pval +"] less than threshold [" + threshold +"] " +site.toStringWithoutGenotypes());
			return true;
		}
		return false;
	}

	/**
	 * Calculates HWE pvalue for the site.  If the site contains no data, return null.
	 * @param site
	 * @return
	 */
	public static Double getHardyWeinbergPvalue (final VariantContext site) {
		ObjectCounter<Integer> c = new ObjectCounter<>();
		for (String s: site.getSampleNames()) {
			Genotype g = site.getGenotype(s);
			if (g!=null && !g.isNoCall()) {
				if (g.isHomRef()) c.increment(0);
				if (g.isHomVar()) c.increment(2);
				if (g.isHet()) c.increment(1);
			}
		}
		if (c.getTotalCount()==0) return null;
		return HardyWeinbergCalculation.hwCalculate(c.getCountForKey(0), c.getCountForKey(1), c.getCountForKey(2));
	}
	@Override
	public void logFilterResults() {
		String msg = String.format("Filter p-value threshold [%f] records pass [%d] records fail [%d] ",this.threshold, this.getRecordsPassed(), this.getRecordsFailed());  
		log.info(msg);										
	}





}
