package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import org.broadinstitute.dropseqrna.utils.StringInterner;
import org.broadinstitute.dropseqrna.vcftools.SampleAssignmentVCFUtils;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;

/**
 * Given a VCF Iterator, gather both the interval for each SNP, as well as the reference and alternate allele, and the average genotype quality of the donors.
 * 
 * The iterator should be pre-filtered to the set of variants that are useful for the analysis.
 * @see SampleAssignmentVCFUtils Use this for some standard filtering methods used in analysis. 
 * @author nemesh
 *
 */
public class SNPInfoCollection {

	private Map<Interval, String> refAllele;
	private Map<Interval, String> altAllele;
	private Map<Interval, Double> averageGQ;
	
	private IntervalList intervalList;
	
	/**
	 * Get the reference allele string
	 * @param snpInterval The snp interval to query
	 * @return The reference allele
	 */
	public String getRefAllele(Interval snpInterval) {		
		return refAllele.get(snpInterval);
	}
	
	/**
	 * Get the alternate allele string
	 * @param snpInterval The snp interval to query
	 * @return The alternatee allele
	 */
	public String getAltAllele(Interval snpInterval) {
		return altAllele.get(snpInterval);
	}

	/**
	 * Get the mean genotype quality of all donors at this SNP site 
	 * @param snpInterval The snp interval to query
	 * @return The mean GQ across donors
	 */
	public Double getAverageGQ(Interval snpInterval) {
		return averageGQ.get(snpInterval);
	}

	/**
	 * Get the list of intervals.
	 * @return The interval list for all sites in the VCF
	 */
	public IntervalList getIntervalList() {
		return intervalList;
	}
	
	/**
	 * Get a map from the SNP Interval to the reference allele for that SNP
	 * @return
	 */
	public Map<Interval, String> getRefAllele() {
		return refAllele;
	}

	/**
	 * Get a map from the SNP Interval to the alternate for that SNP.  Missing alleles are coded as N.
	 * @return
	 */
	public Map<Interval, String> getAltAllele() {
		return altAllele;
	}

	/**
	 * Get a map from the SNP Interval to the average genotype quality for a variant across all included samples.
	 * @return
	 */
	public Map<Interval, Double> getAverageGQ() {
		return averageGQ;
	}
	
	public boolean isEmpty() {
		return this.intervalList.getIntervals().isEmpty();
	}

	/**
	 * For a given VCF iterator, extract the interval list, average genotype quality (GQ), reference and alternate alleles
	 * @param vcfIterator The iterator to extract VCF records from
	 * @param sd The sequence dictionary to use for the interval list
	 * @param preserveIntervalNames Set to true to use the genotype site ID (via site.getID) as the interval name
	 * @param vcfWriter (Optional) write the VCF records from the iterator to this file.  Set to null to ignore.
	 * @param log Write logging to this object if requested.  Set to null to disable logging.
	 */
	public SNPInfoCollection (final Iterator<VariantContext> vcfIterator, final SAMSequenceDictionary sd, boolean preserveIntervalNames, final VariantContextWriter vcfWriter, final Log log) {
		refAllele=new HashMap<>(); 
		altAllele=new HashMap<>();
		averageGQ = new HashMap<>();
		
		SAMFileHeader h = new SAMFileHeader();
		h.setSequenceDictionary(sd);
		IntervalList result = new IntervalList(h);
		if (log!=null) log.info("Scanning VCF to find potential SNP sites");
		
		// cache the reference and alternate alleles
		StringInterner  alleleCache = new StringInterner();
		
		while (vcfIterator.hasNext()) {
			VariantContext site = vcfIterator.next();
			
			Interval variant;
			if (preserveIntervalNames)
				variant = new Interval(site.getContig(), site.getStart(),site.getEnd(), true, site.getID());
			else
				variant = new Interval(site.getContig(), site.getStart(),site.getEnd());
			result.add(variant);
			
			String ref = alleleCache.intern(site.getReference().getBaseString());			
			refAllele.put(variant, ref);
			
			String alt = alleleCache.intern(getAltAllele(site));			
			altAllele.put(variant, alt);
			
			double gqAverage = site.getGenotypes().stream().mapToInt(x -> x.getGQ()).average().orElse(-1);
			averageGQ.put(variant, gqAverage);
			
			// optionally write the VCF record to the writer
			if (vcfWriter!=null)
				vcfWriter.add(site);
			
		}
		if (log!=null) log.info("Found [" + result.getIntervals().size() +"] potential SNP sites to query.");
		this.intervalList=result;
	}
	
	private String getAltAllele (VariantContext site) {
		Allele alt = site.getAltAlleleWithHighestAlleleCount();

		char altBase='N';
		if (alt!=null) {
			byte [] altBases = alt.getBases();
			if (altBases.length>0)
				altBase=StringUtil.byteToChar(altBases[0]);
		}
		return (Character.toString(altBase));
	}
}

