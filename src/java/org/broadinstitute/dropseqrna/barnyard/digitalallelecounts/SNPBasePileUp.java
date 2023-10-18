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
package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.descriptive.moment.Mean;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.StringUtil;

public abstract class SNPBasePileUp implements Locatable,SNPIntervalRecordI {

	private final Interval snpInterval;
	private List<Byte> bases;
	private List<Byte> qualities;

	public SNPBasePileUp (final Interval snpInterval) {
		this.snpInterval=snpInterval;
		this.bases=new ArrayList<>();
		this.qualities=new ArrayList<>();
	}

	public String getSnpID() {
		return this.snpInterval.getName();
	}

	public String getContig() {
		return this.snpInterval.getContig();
	}

	public int getStart() {
		return this.snpInterval.getStart();
	}
	
	public int getEnd () {
		return this.snpInterval.getEnd();
	}

	public Interval getSNPInterval () {
		return this.snpInterval;
	}

	public List<Byte> getBases() {
		return bases;
	}

	public int getNumBases() {
		return this.bases.size();
	}

	public int getCountBase (final char base) {
		Byte baseB = StringUtil.charToByte(base);
		int count=0;
		for (Byte b: bases)
			if (b.equals(baseB)) count++;
		return count;
	}


	public List<Character> getBasesAsCharacters() {
		List<Character> characters = new ArrayList<>(this.bases.size());
		for (Byte b: this.bases)
			characters.add(StringUtil.byteToChar(b));
		return characters;
	}

	public List<Byte> getQualities() {
		return qualities;
	}
	
	public Double getAverageQuality () {
		if (this.qualities.size()==0) return null;
		Mean m = new Mean();
		double [] q = this.qualities.stream().mapToDouble(x -> new Double(x)).toArray();
		return m.evaluate(q);		
	}

	/**
	 * Overwrite any existing bases and qualities with the provided data.
	 *
	 * @param bases
	 * @param qualities
	 */
	public void setBasesAndQualities(final List<Byte> bases, final List<Byte> qualities) {
		if (bases.size()!=qualities.size())
			throw new IllegalArgumentException("Size of bases [" + bases.size() +"] and qualities [" + qualities.size() +"] must be the same.");
		this.bases=bases;
		this.qualities=qualities;
	}

	public void setBasesAndQualities(final byte [] bases, final byte [] qualities) {
		List<Byte> baseList = new ArrayList<>(Arrays.asList(ArrayUtils.toObject(bases)));
		List<Byte> qualList = new ArrayList<>(Arrays.asList(ArrayUtils.toObject(qualities)));
		setBasesAndQualities(baseList, qualList);
	}

	/**
	 * Add a single base and quality score to the pileup.
	 * @param base
	 * @param quality
	 */
	public void addBaseAndQuality (final byte base, final byte quality) {
		this.bases.add(base);
		this.qualities.add(quality);
	}

	/**
	 * For a read on this pileup, get the base and quality of the base that is at the same
	 * position as the SNP.  If the read does not overlap the interval, then return an empty array.
	 * @param r The read
	 * @return A length 2 byte array containing the base and quality.  Empty if the read does not overlap the interval.
	 */
	public byte [] getBaseAndQualityOverlappingInterval (final SAMRecord r) {
		byte [] result = new byte [2];

		List<CigarElement> elements = r.getCigar().getCigarElements();
		Iterator<AlignmentBlock> blocks = r.getAlignmentBlocks().iterator();

		int finalSNPLocalPosition=-1;
		int lengthTraversed=0;

		for (CigarElement ce: elements) {
			CigarOperator co = ce.getOperator();
			// you're in an alignment block
			if (co==CigarOperator.M) {
				// get the next alignment block
				AlignmentBlock b = blocks.next();
				int refStart = b.getReferenceStart();
				int snpLocalPos=this.getStart() - refStart +1;
				int blockLength=b.getLength();

				// is the local position inside this alignment block?
				// if not, move onto the next block.
				if (snpLocalPos >0 && snpLocalPos<=blockLength) {
					// found it!  Done.
					finalSNPLocalPosition=snpLocalPos+lengthTraversed;
					break;
				}
			}
			// consume the bases if necessary and move on to the next element.
			if (co.consumesReadBases())
				lengthTraversed+=ce.getLength();
		}

		// if the position is assigned, then add to the pileup.
		if (finalSNPLocalPosition!=-1) {
			// arrays 0 based.
			byte [] quals = r.getBaseQualities();
			byte [] bases = r.getReadBases();
			byte base = bases[finalSNPLocalPosition-1];
			// char baseChar = StringUtil.byteToChar(base);
			byte qual = quals[finalSNPLocalPosition-1];
			result[0]=base;
			result[1]=qual;
			return (result);
		}
		return (ArrayUtils.EMPTY_BYTE_ARRAY);
	}

	public abstract void addRead (final SAMRecord r);

}
