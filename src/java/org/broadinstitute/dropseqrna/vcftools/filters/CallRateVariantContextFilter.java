/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */
package org.broadinstitute.dropseqrna.vcftools.filters;

import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;

import java.util.*;


/**
 * Filters variants on the number of samples that have a GQ >= genotypeThreshold, or genotypes that have a "no call" set.
 * @author nemesh
 *
 */
public class CallRateVariantContextFilter extends FilteredIterator <VariantContext> {

	private final Log log = Log.getInstance(CallRateVariantContextFilter.class);

	private final int genotypeThreshold;
	private final double fractionPassing;
	private final Set<String> samples;
	private final boolean verbose;

	public CallRateVariantContextFilter (final Iterator<VariantContext> underlyingIterator, final int genotypeThreshold, final double fractionPassing) {
		this (underlyingIterator, genotypeThreshold, fractionPassing, Collections.emptySet(), false);
	}

	public CallRateVariantContextFilter (final Iterator<VariantContext> underlyingIterator, final int genotypeThreshold, final double fractionPassing, final Collection<String> samples) {
		this(underlyingIterator, genotypeThreshold, fractionPassing, samples, false);
	}
	public CallRateVariantContextFilter (final Iterator<VariantContext> underlyingIterator, final int genotypeThreshold, final double fractionPassing, final Collection<String> samples, final boolean verbose) {
		super(underlyingIterator);
		this.genotypeThreshold=genotypeThreshold;
		this.fractionPassing=fractionPassing;
		this.samples=new HashSet<>(samples);
		this.verbose=verbose;
	}


	@Override
	public boolean filterOut(final VariantContext site) {
		Iterator<Genotype> iter = site.getGenotypes().iterator();

		int total=0;
		int pass=0;
		while (iter.hasNext()) {
			Genotype g = iter.next();
			if (this.samples.isEmpty() || this.samples.contains(g.getSampleName())) {
				int quality = g.getGQ();
				if (quality>=this.genotypeThreshold && !g.isNoCall())
					pass++;
				total++;
			}
		}
		double passFrac = (double) pass / (double) total;
		if (passFrac<fractionPassing) {
			if (verbose) log.info("Rejecting variant call rate [" + passFrac +"] less than threshold [" + fractionPassing +"] " +site.toStringWithoutGenotypes());
			return true;
		}
		return false;
	}

	@Override
	public void logFilterResults() {
		String msg = String.format("GQ threshold [%d] fraction donors passing [%f] records pass [%d] records fail [%d] ",this.genotypeThreshold, this.fractionPassing, this.getRecordsPassed(), this.getRecordsFailed());  
		log.info(msg);
		
	}


}
