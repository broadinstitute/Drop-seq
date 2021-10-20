package org.broadinstitute.dropseqrna.utils.readiterators;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.broadinstitute.dropseqrna.utils.CountChangingIteratorWrapper;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.CloseableIterator;
import picard.annotation.LocusFunction;

/**
 * Adds the gene function annotations to reads on the fly. This is useful for
 * DNA reads that need to be evaluated by RNASeq methods.
 * 
 * @author nemesh
 *
 */
public class IgnoreGeneAnnotationTagger extends CountChangingIteratorWrapper<SAMRecord> implements CloseableIterator<SAMRecord> {

	private final String geneTag;
	private final String strandTag;
	private final String functionTag;
	private final Set<LocusFunction> acceptedFunctions;
	private final boolean splitContigByStrand;
	private final String umiTag;
	private final boolean forceRetag;

	/**
	 * Retag the records from the BAM on the fly to convert DNA-Seq data that lacks
	 * gene tags and possibly molecular barcodes. Also adds tags to RNASeq reads
	 * that currently lack tags. Does NOT replace existing tags if the read function
	 * matches the accepted locus functions unless explicitly requested.
	 * 
	 * This adds a gene name tag that matches the chromosome of the read and
	 * (optionally, given splitContigByStrand) the strand of the read, such that the
	 * gene name of the read is set or replaced to be contig/strand. Eg: chr
	 * (splitContigByStrand=false) or chr1+ (splitContigByStrand=true)
	 * 
	 * The gene function is replaced by the first function in the acceptedFunctions
	 * list, such that all reads will pass that filter.
	 * 
	 * If umiTag is not set to null, the read will be tagged with a UMI sequence of
	 * the read name, guaranteeing all read pairs from a BAM come from unique
	 * molecular barcodes.
	 * 
	 * @param underlyingIterator  The BAM file iterator
	 * @param geneTag             The BAM tag used to denote the gene name
	 * @param strandTag           The BAM tag used to denote the gene strand
	 * @param functionTag         The BAM tag used to denote the gene function
	 * @param acceptedFunctions   The list of accepted functions for the data. Reads
	 *                            will be tagged with the first function from the
	 *                            list.
	 * @param splitContigByStrand Should the gene names include the strand of the
	 *                            read?
	 * @param umiTag              If set to a not null value, then the UMI tag will
	 *                            be updated with the read name.
	 * @param forceRetag          Even if the read's locus function matches the
	 *                            acceptedFunctions list, force the read to be
	 *                            retagged.
	 */
	public IgnoreGeneAnnotationTagger(Iterator<SAMRecord> underlyingIterator, String geneTag, String strandTag, String functionTag,
			List<LocusFunction> acceptedFunctions, boolean splitContigByStrand, String umiTag, boolean forceRetag) {
		super(underlyingIterator);
		this.geneTag = geneTag;
		this.strandTag = strandTag;
		this.functionTag = functionTag;
		this.acceptedFunctions = new HashSet<LocusFunction>(acceptedFunctions);
		this.splitContigByStrand = splitContigByStrand;
		this.umiTag = umiTag;
		this.forceRetag=forceRetag;
	}

	/**
	 * Retag the records from the BAM on the fly to convert DNA-Seq data that lacks
	 * gene tags and possibly molecular barcodes. Also adds tags to RNASeq reads
	 * that currently lack tags. Does NOT replace existing tags if the read tag
	 * matches the locus function.
	 * 
	 * This adds a gene name tag that matches the chromosome of the read and
	 * (optionally, given splitContigByStrand) the strand of the read, such that the
	 * gene name of the read is set or replaced to be contig/strand. Eg: chr
	 * (splitContigByStrand=false) or chr1+ (splitContigByStrand=true)
	 * 
	 * @param underlyingIterator  The BAM file iterator
	 * @param geneTag             The BAM tag used to denote the gene name
	 * @param strandTag           The BAM tag used to denote the gene strand
	 * @param functionTag         The BAM tag used to denote the gene function
	 * @param acceptedFunctions   The list of accepted functions for the data. Reads
	 *                            will be tagged with the first function from the
	 *                            list.
	 * @param splitContigByStrand Should the gene names include the strand of the
	 *                            read?
	 * @param forceRetag          Even if the read's locus function matches the
	 *                            acceptedFunctions list, force the read to be
	 *                            retagged.
	 */

	public IgnoreGeneAnnotationTagger(Iterator<SAMRecord> underlyingIterator, String geneTag, String strandTag, String functionTag,
			List<LocusFunction> acceptedFunctions, boolean splitContigByStrand, boolean forceRetag) {
		this(underlyingIterator, geneTag, strandTag, functionTag, acceptedFunctions, splitContigByStrand, null, forceRetag);
	}

	@Override
	/**
	 * Adds new tags (or over writes old ones) such that each read now has a gene
	 * name for the contig the reads is on, with the strand of the read appended.
	 * 
	 */
	protected void processRecord(SAMRecord rec) {
		// If the UMI tag is provided, use the read name as the UMI.
		if (umiTag != null)
			rec.setAttribute(this.umiTag, rec.getReadName());

		String functionList = rec.getStringAttribute(this.functionTag);
		// if locus function is null, there's no tag for function on the read and we can
		// add tags.
		// Otherwise, check if the tag on the read is an "accepted" tag. If so, keep it.
		if (functionList != null) {
			final LocusFunction[] locusFunctions = (functionList == null ? null : GeneFunctionIteratorWrapper.getLocusFunctionFromRead(functionList));
			Set<LocusFunction> locusSet = new HashSet<>(Arrays.asList(locusFunctions));
			// if the read has locus functions that were requested for use, then rely on
			// them.
			locusSet.retainAll(acceptedFunctions);
			if (locusSet.size() > 0) {
				queueRecordForOutput(rec);
				return;
			}
		}

		// We do encode the strand in the contig name
		// Disabled for now because it can double count UMIs.
		if (splitContigByStrand) {
			// Annotate the contig by the strand, effectively splitting the contig into 2
			// "genes"
			if (rec.getReadNegativeStrandFlag()) {
				rec.setAttribute(geneTag, rec.getContig() + "-");
			} else {
				rec.setAttribute(geneTag, rec.getContig() + "+");
			}
		} else {
			// leave the contig as a single "gene"
			rec.setAttribute(geneTag, rec.getContig());
		}

		// Set the strand to be the same as the read.
		if (rec.getReadNegativeStrandFlag())
			rec.setAttribute(strandTag, "-");
		else
			rec.setAttribute(strandTag, "+");

		// set the record as having an accepted tag.
		String thisFunctionTag = acceptedFunctions.iterator().next().toString();
		// TODO: at some point couple this more tightly to the messaging of the main
		// class. If this implementation changes, must change that as well.
		rec.setAttribute(functionTag, thisFunctionTag);
		queueRecordForOutput(rec);
	}

	public String getGeneTag() {
		return geneTag;
	}

	public String getStrandTag() {
		return strandTag;
	}

	public String getFunctionTag() {
		return functionTag;
	}

	public String getUmiTag() {
		return umiTag;
	}

}
