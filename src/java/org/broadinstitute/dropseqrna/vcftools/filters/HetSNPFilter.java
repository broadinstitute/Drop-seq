package org.broadinstitute.dropseqrna.vcftools.filters;

import java.util.Iterator;

import org.broadinstitute.dropseqrna.utils.FilteredIterator;

import htsjdk.samtools.util.Log;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;

public class HetSNPFilter extends FilteredIterator <VariantContext>  {
	
	private static final Log log = Log.getInstance(HetSNPFilter.class);
	
	public HetSNPFilter(final Iterator<VariantContext> underlyingIterator) {
		super(underlyingIterator);
	}

	@Override
	public boolean filterOut(final VariantContext rec) {
		for (Genotype g : rec.getGenotypes()) {
			if (!g.isHet()) 
				return true;
		}
		return false;
	}

	@Override
	public void logFilterResults() {
		String msg = String.format("Records pass [%d] records fail [%d] ",this.getRecordsPassed(), this.getRecordsFailed());  
		log.info(msg);										
	}

}

