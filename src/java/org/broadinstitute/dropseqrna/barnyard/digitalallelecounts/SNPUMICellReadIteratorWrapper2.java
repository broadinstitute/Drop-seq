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
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.broadinstitute.dropseqrna.barnyard.Utils;
import org.broadinstitute.dropseqrna.utils.CountChangingIteratorWrapper;
import org.broadinstitute.dropseqrna.utils.IntervalTagComparator;
import org.broadinstitute.dropseqrna.utils.PassFailTrackingIteratorI;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.OverlapDetector;

public class SNPUMICellReadIteratorWrapper2 extends CountChangingIteratorWrapper<List<SAMRecord>> {
	private static final Log log = Log.getInstance(SNPUMICellReadIteratorWrapper.class);
	
	private String cellBarcodeTag;
	private Set<String> cellBarcodeList;
	private String geneTag;
	private final String snpTag;
	private final OverlapDetector<Interval> snpIntervals;
	private final Map<Interval, Double> meanGenotypeQuality;
	
	private int recordsPass=0;
	private int recordsTotal=0;
	private boolean isLogged=false;
	private int failFastThreshold=-1;
	
	
	/**
	 * Filters and copies reads for generating SNPCellUMIBasePileups.
	 * Filters out reads in which the read cell barcode does not match a barcode in the list (if list is not null)
	 * Reads that are marked as secondary or supplementary are filtered out
	 * Filters reads based on read map quality, removing reads below that quality
	 * Optionally filters reads where the annotated gene and the strand of the read don't match, or can clone a read and return it multiple times
	 * if the read maps to more than one gene and <assignReadsToAllGenes> is true.
	 * @param underlyingIterator Source of reads.  These reads are sorted by cell/gene/umi then grouped, such that each List contains all the reads for a "UMI".
	 * @param cellBarcodeTag The cell barcode BAM tag
	 * @param cellBarcodeList A list of cell barcodes, or null to ignore.  If populated, reads where the cellBarcodeTag matches one of these Strings will be retained
	 * @param geneTag The gene/exon tag.
	 * @param strandTag The strand tag
	 * @param readMQ The minimum map quality of a read to be retained.
	 * @param A map from SNP Intervals to the mean genotype quality.  In cases where a read has 
	 * multiple SNPs, the SNP with the highest average genotype quality will be selected to avoid read "double counting".  
	 */
	public SNPUMICellReadIteratorWrapper2(final Iterator<List<SAMRecord>> underlyingIterator,
                                         final IntervalList snpIntervals,
                                         final String cellBarcodeTag,
                                         final Collection<String> cellBarcodeList,
                                         final String geneTag,
                                         final String snpTag,
                                         final int readMQ,
                                         final Map<Interval, Double> meanGenotypeQuality) {
        
		super(underlyingIterator);
		
		if (meanGenotypeQuality==null)
			log.warn("No Mean Genotype Quality supplied, will not filter multiple variants on the same UMI");
		
		this.cellBarcodeTag = cellBarcodeTag;
		this.cellBarcodeList = new HashSet<String>(cellBarcodeList);
		this.geneTag=geneTag;
		this.snpTag=snpTag;
		this.meanGenotypeQuality=meanGenotypeQuality;

		// construct OverlapDetector
		OverlapDetector<Interval> od = new OverlapDetector<>(0, 0);
		od.addAll(snpIntervals.getIntervals(), snpIntervals.getIntervals());
		this.snpIntervals=od;
	}
	
	public void setFailFastThreshold (int count) {
		this.failFastThreshold=count;
	}
	
	

    @Override
    protected void processRecord(final List<SAMRecord> recList) {    	
    	SAMRecord r=recList.get(0);
        String cellBC=r.getStringAttribute(cellBarcodeTag);
        String geneList = r.getStringAttribute(this.geneTag);        

        // if there are cell barcodes to filter on, and this read's cell barcode isn't one of them, then move on to the next read;
        if (this.cellBarcodeList!=null && !cellBarcodeList.contains(cellBC))
			return;
        // if the read have any genes.
        if (geneList==null)
			return;
        processGene(recList);
    }

	/**
	 * Check if any reads of a UMI (cell/gene/molecular barcode) overlap any SNPs in the OverlapDetector.  Tag reads with SNPs.
	 * If more than 1 SNP tags a UMI, select the highest quality SNP for that UMI.
	 * Simplified since data goes through GeneFunctionIteratorWrapper to take care of how reads/genes interact.
	 */
	private void processSNP (final List<SAMRecord> recList) {
		checkFastFail();
		
		recordsTotal++;
		
		Map<Interval, Set<SAMRecord>> snpIntervalToReadMap = new HashMap<>();
		// find all 
		for (SAMRecord r : recList) {
			List<AlignmentBlock> blocks = r.getAlignmentBlocks();
			for (AlignmentBlock b: blocks) {
				int start = b.getReferenceStart();
				int end = start + b.getLength() -1;

				Interval i = new Interval(r.getReferenceName(), start, end);			
				Collection<Interval> overlaps = this.snpIntervals.getOverlaps(i);
				for (Interval interval: overlaps) {
					Set<SAMRecord> recs = snpIntervalToReadMap.get(interval);
					if (recs==null) {
						recs = new HashSet<>();						
					}
					recs.add(r);
					snpIntervalToReadMap.put(interval, recs);
				}
			}
		}	

		// exit early if no SNPs found.
		if (snpIntervalToReadMap.size()==0) {
			return;
		}
								
		// exit without selecting a SNP if there's only a single SNP.
		if (snpIntervalToReadMap.size()==1) {
			String intervalString=IntervalTagComparator.toString(snpIntervalToReadMap.keySet().iterator().next());
			for (SAMRecord r: recList) {				
				r.setAttribute(this.snpTag, intervalString);					
			}
			queueRecordForOutput(recList);
			this.recordsPass++;
			return;
		}
		
		// if there are multiple SNPs seen by this UMI, select the best SNP, and the reads that support that SNP.
		if (this.meanGenotypeQuality!=null) {					
			Interval snpInterval = getBestSNP(snpIntervalToReadMap.keySet());
			String intervalString=IntervalTagComparator.toString(snpInterval);
			
			List<SAMRecord> result = new ArrayList<> (snpIntervalToReadMap.get(snpInterval));			
			for (SAMRecord r: result) {
				r.setAttribute(this.snpTag, intervalString);	
			}			
			queueRecordForOutput(result);
			this.recordsPass++;
			return;
		}
		
		// Return one copy of each read for each SNP that touches it.
		List<SAMRecord> result = new ArrayList<>(); 
		for (Interval snp:snpIntervalToReadMap.keySet()) {
			Collection<SAMRecord> recs = snpIntervalToReadMap.get(snp);
			for (SAMRecord r: recs) {
				SAMRecord rr = Utils.getClone(r);
				rr.setAttribute(this.snpTag, IntervalTagComparator.toString(snp));
				result.add(rr);
			}									
		}
		queueRecordForOutput(result);						
		this.recordsPass++;
	}
	
	private Interval getBestSNP (Collection<Interval> snpIntervals) {
		// set the worse score to be worse than the missing value of -1.
		double maxGQ=-2d;
		Interval best=null;
		
		for (Interval i: snpIntervals) {
			double gq = this.meanGenotypeQuality.get(i);
			if (gq>maxGQ) {
				maxGQ=gq;
				best=i;
			}
		}
		return (best);
	}


	/**
	 * For a read, check and see if the read maps to more than 1 gene, and if the gene matches the strand (optional, controlled by useStrandInfo flag).
	 * If useStrandInfo is true and the read strand does not match the gene strand, discard the read.
	 * If useStrandInfo is true and the read strand matches the gene strand, keep the read.
	 * If useStrandInfo is false, keep the read
	 * If there is more than 1 gene tagging the read, make a copy of the read for each gene.  This also uses the strand option.
	 */
	private void processGene (final List<SAMRecord> recList) {
		// split all reads into individual genes.
		List<SAMRecord> result = new ArrayList<>();
		
		for (SAMRecord r: recList) {
			String geneList = r.getStringAttribute(this.geneTag);
			String [] genes = geneList.split(",");
			// if there's one gene, no need to split
			if (genes.length==1) {
				result.add(r);
			} else { // otherwise make clones of reads for each gene, usually the case with ambiguous genes like intron/intron gene pairs.
				for (String g : genes) {
					SAMRecord rr = Utils.getClone(r);
					rr.setAttribute(geneTag, g);
					result.add(rr);
				}
			}
			
		}
		processSNP(result);
		
	}

		
	@Override
    public boolean hasNext() {
		// checkFastFail();
        boolean hasNext = super.hasNext();
        if (!hasNext & !isLogged) {
        	String msg = String.format("UMIs that with at least one SNP [%d] of total UMIs [%d] ",this.recordsPass, this.recordsTotal);  
    		log.info(msg);    		
    		//TODO: why does hasNext get called after it returns null?
    		isLogged=true;
        }
        return hasNext;
    }
	
	private void checkFastFail () {
		if (this.failFastThreshold!=-1 && this.recordsTotal>= this.failFastThreshold & this.recordsPass==0) {
			String msg = String.format("Did not encounter any transcribed SNPs after testing [%d] UMIs", this.failFastThreshold);   
			log.error(msg);
			throw new IllegalStateException();
		}
	}

}
