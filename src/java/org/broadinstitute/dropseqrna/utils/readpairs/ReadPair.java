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
package org.broadinstitute.dropseqrna.utils.readpairs;

import htsjdk.samtools.SAMRecord;

import org.broadinstitute.dropseqrna.junctionlibrary.JunctionSamUtils;

public class ReadPair {

		private SAMRecord read1=null;
		private SAMRecord read2=null;
		
		/**
		public ReadPair (SAMRecord read1) {
			this.read1=read1;
		}
		*/
		
		public ReadPair (SAMRecord read1, SAMRecord read2) {
			if (read1.getFirstOfPairFlag()) {
				this.read1 = read1;
			}
			if (read1.getSecondOfPairFlag()) {
				this.read2=read1;
			}
			if (read2.getFirstOfPairFlag()) {
				this.read1=read2;
			}
			if (read2.getSecondOfPairFlag()) {
				this.read2=read2;
			}
		}
		
		
		
		
		public ReadPair getFlippedReadPair () {
			ReadPair p = new ReadPair();
			p.setRead2(this.getRead1());
			p.setRead1(this.getRead2());
			return (p);
		}
		
		public ReadPair() {}
		
		public SAMRecord getRead1() {
			return read1;
		}

		public void setRead1(SAMRecord read1) {
			this.read1 = read1;
		}

		public SAMRecord getRead2() {
			return read2;
		}

		public void setRead2(SAMRecord read2) {
			this.read2 = read2;
		}
		
		public SAMRecord getFirstRead() {
			if (read1.getFirstOfPairFlag()) {
				return this.read1;
			}
			if (read2.getFirstOfPairFlag()) {
				return this.read2;
			}
			return null;
		}
		
		public SAMRecord getSecondRead() {
			if (read1.getFirstOfPairFlag()) {
				return this.read2;
			}
			if (read2.getFirstOfPairFlag()) {
				return this.read1;
			}
			return null;
		}
		
		public SAMRecord getLeftRead () {
			if (this.read1.getAlignmentStart()<=this.read2.getAlignmentStart()) {
				return this.read1;
			}
			return this.read2;
		}
		
		public SAMRecord getRightRead () {
			if (this.read2.getAlignmentStart()>=this.read1.getAlignmentStart()) {
				return this.read2;
			}
			return this.read1;
		}
		
		public boolean testProperlyPaired() {
			 return (JunctionSamUtils.getInstance().testPairedRead(this.read1, this.read2));
		}
		
				
		public String toString () {
			StringBuilder b = new StringBuilder();
			b.append(samToString(this.read1)+"\n");
			b.append(samToString(this.read2)+"\n");
			return b.toString();
		}
		
		private String samToString (SAMRecord r) {
			StringBuilder b = new StringBuilder();
			b.append(" read name " + r.getReadName());
			b.append(" is first read " + r.getFirstOfPairFlag());
			b.append(" is unmapped " + r.getReadUnmappedFlag());
			b.append(" map quality " +r.getMappingQuality());
			b.append(" reference " + r.getReferenceName());
			b.append(" map start pos " + r.getAlignmentStart());
			b.append(" cigar " + r.getCigarString());
			b.append(" mate unmapped " + r.getMateUnmappedFlag());
			b.append(" mate reference " + r.getMateReferenceName());
			b.append(" mate start pos " + r.getMateAlignmentStart());
			return b.toString();	
		}
}
