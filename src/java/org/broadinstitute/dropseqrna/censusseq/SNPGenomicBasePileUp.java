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

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Interval;

import java.util.Collection;
import java.util.List;

import org.broadinstitute.dropseqrna.barnyard.digitalallelecounts.SNPBasePileUp;

public class SNPGenomicBasePileUp extends SNPBasePileUp {

	private boolean finalized;

	private ListMultimap<String, SAMRecord> recs = ArrayListMultimap.create();

	public SNPGenomicBasePileUp (final Interval snpInterval) {
		super(snpInterval);
		this.finalized=false;
	}

	/**
	 * Create a new pileup from this object where bases in the pileup are filtered by a minimum quality threshold.
	 * @param minQuality
	 * @return
	 */
	public SNPGenomicBasePileUp getFilteredPileupByBaseQuality (final byte minQuality) {
		if (!this.finalized) finalizePileup();
		SNPGenomicBasePileUp r = new SNPGenomicBasePileUp(this.getSNPInterval());
		List<Byte> bases = getBases();
		List<Byte> quals = getQualities();
		for (int i=0; i<bases.size(); i++) {
			Byte qual = quals.get(i);
			if (qual>=minQuality)
				r.addBaseAndQuality(bases.get(i), qual);
		}
		r.finalized=true;
		return r;
	}



	@Override
	/**
	 * Adds a read to the pileup.
	 * Paired end genomic reads are a little more complicated than single ended reads.
	 * We hold all reads in memory until we declare we're done adding reads, then find read pairs and
	 * add data a pair at a time.
	 * If 2 reads in a pair overlap a SNP, we check for consistency - if both reads agree on the allele,
	 * that allele is added to the pileup.  If the reads of the pair disagree at that base, we discard the pair.
	 *
	 * If this pileup has been finalized, no more reads can be added.
	 * @param r the read to add to the pileup.
	 */
	public void addRead(final SAMRecord r) {
		if (this.finalized) throw new IllegalStateException("Adding a read to an already finalized pileup!");
		recs.put(r.getReadName(), r);
	}

	public void addReads(final Collection<SAMRecord> recs) {
		for (SAMRecord r: recs)
			addRead(r);
	}

	/**
	 * Add all the reads that have previously been added to the final pileup.
	 * After this call, no more reads can be added to this pileup object.
	 *
	 */
	public void finalizePileup () {
		this.finalized=true;
		for (String recName : recs.keySet()) {
			List<SAMRecord> sr = recs.get(recName);
			addRecordPairToPileUp(sr);
		}
	}

	/**
	 * Adds all reads with a read name to the pileup.
	 * This should be one or two reads.
	 * If one read, simply accept the read if it overlaps the interval
	 * If two reads, only add if the base from both reads agrees, and if so take the higher base quality.
	 * @param recs
	 */
	void addRecordPairToPileUp (final List<SAMRecord> recs) {
		// the easy case of a single read.
		if (recs.size()==1) {
			byte [] baseAndQual = getBaseAndQualityOverlappingInterval(recs.get(0));
			if (baseAndQual.length!=0)
				addBaseAndQuality(baseAndQual[0], baseAndQual[1]);
			return; // short circuit out.
		}
		// paired read.
		if (recs.size()==2) {
			byte [] baseAndQualOne = getBaseAndQualityOverlappingInterval(recs.get(0));
			byte [] baseAndQualTwo = getBaseAndQualityOverlappingInterval(recs.get(1));
			if (baseAndQualOne[0]==baseAndQualTwo[0]) {
				if (baseAndQualOne[1] > baseAndQualTwo[1])
					addBaseAndQuality(baseAndQualOne[0], baseAndQualOne[1]);
				else
					addBaseAndQuality(baseAndQualOne[0], baseAndQualTwo[1]);
			} else
				return; // exit because the bases called on two reads overlapping the same SNP are different.
		} else // whaaaat?
			throw new IllegalStateException("Should never have more than 2 reads in a read pair.");
	}

	@Override
	public String toString () {
		StringBuilder b = new StringBuilder();
		b.append("snp [" + this.getSNPInterval().toString() +"] ");
		b.append(this.getBasesAsCharacters().toString()+ " ");
		b.append(this.getQualities().toString());
		return b.toString();
	}



}
