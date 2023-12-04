package org.broadinstitute.dropseqrna.barnyard.digitalallelecounts;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.StringInterner;
import org.broadinstitute.dropseqrna.vcftools.SampleAssignmentVCFUtils;

import java.util.*;

import static org.broadinstitute.dropseqrna.vcftools.filters.MinorAlleleFreqVariantContextFilter.calculateMinorAlleleFrequency;


/**
 * Given a VCF Iterator, gather both the interval for each SNP, as well as the reference and alternate allele, and the average genotype quality of the donors.
 * <p>
 * The iterator should be pre-filtered to the set of variants that are useful for the analysis.
 * @see SampleAssignmentVCFUtils Use this for some standard filtering methods used in analysis. 
 * @author nemesh
 *
 */
public class SNPInfoCollection {

	private final Map<Interval, String> refAllele;
	private final Map<Interval, String> altAllele;
	private final Map<Interval, Double> averageGQ;

	private final IntervalList intervalList;

	/**
	 * Get the reference allele string
	 *
	 * @param snpInterval The snp interval to query
	 * @return The reference allele
	 */
	public String getRefAllele(Interval snpInterval) {
		return refAllele.get(snpInterval);
	}

	/**
	 * Get the alternate allele string
	 *
	 * @param snpInterval The snp interval to query
	 * @return The alternate allele
	 */
	public String getAltAllele(Interval snpInterval) {
		return altAllele.get(snpInterval);
	}

	/**
	 * Get the mean genotype quality of all donors at this SNP site
	 *
	 * @param snpInterval The snp interval to query
	 * @return The mean GQ across donors
	 */
	public Double getAverageGQ(Interval snpInterval) {
		return averageGQ.get(snpInterval);
	}

	/**
	 * Get the list of intervals.
	 *
	 * @return The interval list for all sites in the VCF
	 */
	public IntervalList getIntervalList() {
		return intervalList;
	}

	/**
	 * Get a map from the SNP Interval to the reference allele for that SNP
	 *
	 * @return Map from interval to the reference allele
	 */
	public Map<Interval, String> getRefAllele() {
		return refAllele;
	}

	/**
	 * Get a map from the SNP Interval to the alternate for that SNP.  Missing alleles are coded as N.
	 *
	 * @return Map from interval to most common alternate allele
	 */
	public Map<Interval, String> getAltAllele() {
		return altAllele;
	}

	/**
	 * Get a map from the SNP Interval to the average genotype quality for a variant across all included samples.
	 *
	 * @return Map from interval to genotype quality
	 */
	public Map<Interval, Double> getAverageGQ() {
		return averageGQ;
	}

	public boolean isEmpty() {
		return this.intervalList.getIntervals().isEmpty();
	}


	/**
	 * For a given VCF iterator, extract the interval list, average genotype quality (GQ), reference and alternate alleles
	 *
	 * @param vcfIterator           The iterator to extract VCF records from
	 * @param sd                    The sequence dictionary to use for the interval list
	 * @param preserveIntervalNames Set to true to use the genotype site ID (via site.getID) as the interval name
	 * @param vcfWriter             (Optional) write the VCF records from the iterator to this file.  Set to null to ignore.
	 * @param log                   Write logging to this object if requested.  Set to null to disable logging.
	 */
	public SNPInfoCollection(final Iterator<VariantContext> vcfIterator, final SAMSequenceDictionary sd, boolean preserveIntervalNames, final VariantContextWriter vcfWriter, final Log log) {
		refAllele = new HashMap<>();
		altAllele = new HashMap<>();
		averageGQ = new HashMap<>();

		SAMFileHeader h = new SAMFileHeader();
		h.setSequenceDictionary(sd);
		IntervalList result = new IntervalList(h);
		if (log != null) log.info("Scanning VCF to find potential SNP sites");

		// cache the reference and alternate alleles
		StringInterner alleleCache = new StringInterner();

		//  Group VariantContext sites by chromosome and position, and select the most representative
		//  alternate allele of the multi-allelic record to represent the pileup alt allele annotation.
		GroupingIterator<VariantContext> groupingIterator = new GroupingIterator<>(vcfIterator, new VariantContextPositionComparator(sd));

		while (groupingIterator.hasNext()) {
			List<VariantContext> listVC = groupingIterator.next();
			// figure out which VC is the representative one.

			VariantContext site = getBestVC(listVC);

			Interval variant;
			if (preserveIntervalNames)
				variant = new Interval(site.getContig(), site.getStart(), site.getEnd(), true, site.getID());
			else
				variant = new Interval(site.getContig(), site.getStart(), site.getEnd());
			result.add(variant);

			String ref = alleleCache.intern(site.getReference().getBaseString());
			refAllele.put(variant, ref);

			String alt = alleleCache.intern(getAltAllele(site));
			altAllele.put(variant, alt);

			double gqAverage = site.getGenotypes().stream().mapToInt(Genotype::getGQ).average().orElse(-1);
			averageGQ.put(variant, gqAverage);

			// optionally write the VCF record to the writer
			if (vcfWriter != null)
				vcfWriter.add(site);

		}
		if (log != null) log.info("Found [" + result.getIntervals().size() + "] potential SNP sites to query.");
		this.intervalList = result;
	}


	private String getAltAllele(VariantContext site) {
		Allele alt = site.getAltAlleleWithHighestAlleleCount();
		char altBase = 'N';
		if (alt != null) {
			byte[] altBases = alt.getBases();
			if (altBases.length > 0)
				altBase = StringUtil.byteToChar(altBases[0]);
		}
		return (Character.toString(altBase));
	}

	/**
	 * Given a list of VariantContext objects that are all at the same genomic position,
	 * pick a representative one for alt allele capture.
	 * This would be the one with the highest minor allele frequency in the population.
	 * @param list
	 * @return
	 */
	private VariantContext getBestVC(List<VariantContext> list) {
		if (list.size()==1)
			return list.get(0);

		List<Double> freqList = new ArrayList<>();
		for (VariantContext vc :  list) {
			double maf = calculateMinorAlleleFrequency(vc, -1);
			maf = Math.min(maf, 1-maf);
			freqList.add(maf);
		}

		int maxIndex = freqList.indexOf(Collections.max(freqList));
		return list.get(maxIndex);
	}
}

class VariantContextPositionComparator implements Comparator<VariantContext> {
	private final SAMSequenceDictionary dict;
	public VariantContextPositionComparator (SAMSequenceDictionary dict) {
		this.dict=dict;
	}
	@Override
	public int compare(final VariantContext o1, final VariantContext o2) {
		int seqIdx1 = dict.getSequenceIndex(o1.getContig());
		int seqIdx2 = dict.getSequenceIndex(o2.getContig());
		int cmp = seqIdx1 - seqIdx2;
		if (cmp==0)
			cmp = o1.getStart() - o2.getStart();
		if (cmp==0)
			cmp = o1.getEnd() - o2.getEnd();
		return cmp;
	}
}

