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
package org.broadinstitute.dropseqrna.censusseq;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.variant.variantcontext.VariantContext;

import org.broadinstitute.dropseqrna.censusseq.JointIteratorCounter;

public class VCFPileupJointIterator {

	private static final Log log = Log.getInstance(VCFPileupJointIterator.class);

	private final PeekableIterator<SNPGenomicBasePileUp> peekablePileUpIter;
	private final PeekableIterator<VariantContext> vcfIterator;
	private final JointIteratorCounter counter;

	private final SAMSequenceDictionary sd;
	private JointResult next=null;

	public VCFPileupJointIterator (final PeekableIterator<SNPGenomicBasePileUp> peekablePileUpIter,
			final PeekableIterator<VariantContext> vcfIterator, final SAMSequenceDictionary sd) {
		this.peekablePileUpIter=peekablePileUpIter;
		this.vcfIterator=vcfIterator;
		this.sd=sd;
		counter = new JointIteratorCounter();
		// prime the iterator.
		getNextSet();
	}

	public boolean hasNext() {
		if (this.next==null) getNextSet(); // iterates until you have a result, or you're out of results.
		return this.next!=null;
	}

	/**
	 * Get the next result.
	 * @return
	 */
	public JointResult next() {
		if (this.next==null) getNextSet(); // iterates until you have a result, or you're out of results.
		// the result you'll return
		JointResult result = this.next;
		// since you are handing out a result, the cached result is null.
		this.next=null;
		return result;
	}



	private void getNextSet() {
		// if you're out of records, report how many total records were processed.
		// if (!vcfIterator.hasNext()|| ! peekablePileUpIter.hasNext())
		// 	log.info("Processed [" + counter.BOTH + "] SNPs in BAM + VCF");
		
		while (vcfIterator.hasNext() && peekablePileUpIter.hasNext()) {
			SNPGenomicBasePileUp pileUp = peekablePileUpIter.peek();
			// only have to compare the first record.
			int cmp = CensusSeqUtils.compareRecords(pileUp, vcfIterator.peek(), sd);

			if (cmp < 0) {
				peekablePileUpIter.next();
				counter.SAMPLE_GENE_ITER++;
				// for debugging...we expect that any SNP pileup should have come with a VCF record...
				cmp = CensusSeqUtils.compareRecords(pileUp, vcfIterator.peek(), sd);
			} else if (cmp > 0) {
				vcfIterator.next();
				counter.VCF_ITER++;
			} else if (cmp == 0) {
				counter.BOTH++;
				counter.SAMPLE_GENE_ITER++;
				counter.VCF_ITER++;
				if (counter.BOTH % 250000 == 0)
					log.info("Processed [" + counter.BOTH + "] SNPs in BAM + VCF");

				// grab the next record and process it.
				VariantContext vc = vcfIterator.next();
				// write the VCF record out if needed.

				pileUp = peekablePileUpIter.next();

				// pileups can contain 0 reads because of filtering of overlapping intervals. Le sigh.
				if (pileUp.getNumBases() > 0) {
					JointResult jr = new JointResult(pileUp, vc);
					this.next=jr;
					break; // you found a result, break the loop.
				} else
					// if there aren't bases for the pileup, this position was rejected.  Pretty rare.
					counter.REJECTED++;
			}
		}		
		
	}

	public JointIteratorCounter getCounter() {
		return counter;
	}

	public void close() {
		this.peekablePileUpIter.close();
		this.vcfIterator.close();
	}

	public class JointResult {

		private final SNPGenomicBasePileUp pileup;
		private final VariantContext vc;

		public JointResult (final SNPGenomicBasePileUp pileup, final VariantContext vc) {
			if ((pileup.getSNPInterval().getStart()!=vc.getStart()) || !pileup.getChromosome().equals(vc.getContig()))
				log.warn("Pileup and VCF at different positions");
			this.pileup=pileup;
			this.vc=vc;
		}

		public SNPGenomicBasePileUp getPileup() {
			return pileup;
		}

		public VariantContext getVc() {
			return vc;
		}





	}



}
