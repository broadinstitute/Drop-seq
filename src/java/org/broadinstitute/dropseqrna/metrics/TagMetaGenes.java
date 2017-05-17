package org.broadinstitute.dropseqrna.metrics;

import java.io.File;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.IntStream;

import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang.builder.CompareToBuilder;
import org.broadinstitute.dropseqrna.barnyard.Utils;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.CustomBAMIterators;
import org.broadinstitute.dropseqrna.utils.FilteredIterator;
import org.broadinstitute.dropseqrna.utils.GroupingIterator;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.ReadNameComparator;
import org.broadinstitute.dropseqrna.utils.SamHeaderUtil;
import org.broadinstitute.dropseqrna.utils.StringInterner;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;
import org.broadinstitute.dropseqrna.utils.readiterators.MapQualityFilteredIterator;
import org.broadinstitute.dropseqrna.utils.readiterators.SamRecordSortingIteratorFactory;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileHeader.SortOrder;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordQueryNameComparator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.StringUtil;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.MathUtil;

@CommandLineProgramProperties(
        usage = "A special case tagger.  Using primary and secondary alignments, try to discover genes that have high homology, where reads map to both genes consistently.",
        usageShort = "Discovery and tag meta-genes",
        programGroup = DropSeq.class
)
public class TagMetaGenes extends CommandLineProgram{

	private final Log log = Log.getInstance(TagMetaGenes.class);
	private ProgressLogger pl = new ProgressLogger(log);

	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze")
	public File INPUT;

	@Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output BAM, written with new Gene/Exon tag")
	public File OUTPUT;

	@Option (doc="A report on which genes are associated to each other.")
	public File REPORT;

	@Option(doc = "The tag to combine.")
	public String GENE_EXON_TAG="GE";

	public String STRAND_TAG="GS";

	/*
	@Option (doc="The metagene tag to add to records")
	public String METAGENE_TAG="MG";
	*/

	@Option(doc="The minimum map quality of the read to get a metatag.")
	public Integer MIN_READ_MQ=3;

	@Option(doc="The maximum map quality of the read to get a metatag.")
	public Integer MAX_READ_MQ=3;

	// @Option(doc="Only consider reads that are within <MAXIMUM_READ_EDIT_DISTANCE> of the place they are mapped to on the reference. ", optional=true)
	private Integer MAXIMUM_READ_EDIT_DISTANCE=null;

	@Option(doc="The map quality of uniquely mapped reads.")
	public Integer UNIQUE_READ_MQ=255;

	@Option (doc="Should reads that don't map to a gene be filtered out of the output report?  These are reads that"
			+ "map to both a gene and a non-gene segment.  Turn this filter on to determine how many multi-mapping reads map"
			+ "both to a gene's exon and a non-exonic segment anywhere else in the genome.")
	public Boolean FILTER_MISSING_GENES=true;

	private String MISSING_GENE="NONE";

	private String DELIMITER=":";

	@Override
	protected int doWork() {
		IOUtil.assertFileIsReadable(this.INPUT);
		IOUtil.assertFileIsWritable(this.REPORT);

		SamReader reader = SamReaderFactory.makeDefault().open(INPUT);

		// time: 92.69 minutes
		// GroupingIterator<SAMRecord> iter = getGroupedReadNameIterator(reader);

		// time: 77.70 minutes
		GroupingIterator<SAMRecord> iter = getGroupedReadNameIteratorFaster(reader);

		SAMFileHeader header = reader.getFileHeader();
		SamHeaderUtil.addPgRecord(header, this);
		SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT);

		ObjectCounter<String> geneClusterCounter = new ObjectCounter<>();
		ObjectCounter<String> uniqueGeneCounter = new ObjectCounter<>();
		ProgressLogger pl = new ProgressLogger(this.log);

		while (iter.hasNext()) {
			List<SAMRecord> records = iter.next();
			for (SAMRecord r: records) pl.record(r);
			geneClusterCounter= populateMultiMapperGenes(records, geneClusterCounter);
			uniqueGeneCounter = populateUniqueMapperGenes(records, uniqueGeneCounter);
		}

		List<MetaGene> metaGenes = buildMetaGenes(geneClusterCounter, uniqueGeneCounter);
		writeReport (metaGenes, this.REPORT, this.FILTER_MISSING_GENES);


		CloserUtil.close(reader);
		// writeReport(geneClusterCounter, uniqueGeneCounter, this.REPORT);

		// here's where you'd write out the newly metagene tagged BAM.

		writer.close();
		return 0;
	}

	private void writeReport (final List<MetaGene> metaGenes, final File outFile, final boolean filterNonGenes) {
		DecimalFormat percentageFormat = new DecimalFormat("###.#");
		PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(outFile));
		// write header
		String [] line={"META_GENE", "META_GENE_COUNT", "UNIQUE_GENE_COUNT", "TOTAL_UNIQUE_READS", "PCT_MULTIMAPPED_READS"};
		String h = StringUtils.join(line, "\t");
        out.println(h);
        for (MetaGene mg: metaGenes) {
        	if (this.FILTER_MISSING_GENES && mg.hasReadsOnNonGenicRegions()) continue; // skip meta-genes that include non-genic regions
        	String uniqueGeneCountsString = StringUtils.join(getString(mg.getUniqueGeneCounts(filterNonGenes)), this.DELIMITER);
        	String [] b ={mg.metaGeneName, Integer.toString(mg.metaGeneReadCount), uniqueGeneCountsString, Integer.toString(mg.getTotalUniqueReads(filterNonGenes)), percentageFormat.format(mg.getPercentMultiMappedReads(filterNonGenes))};
        	String r = StringUtils.join(b, "\t");
            out.println(r);
        }
        out.close();
	}

	private String [] getString (final int [] values) {
		String [] r = new String [values.length];
		for (int i=0; i<values.length; i++)
			r[i]=Integer.toString(values[i]);
		return (r);
	}

	private int [] getUniqueGeneCounts (final String concatonatedGeneList, final String delimiter, final ObjectCounter<String> uniqueGeneCounter) {
		// get the individual gene counts.
    	String [] genes = StringUtils.split(concatonatedGeneList, delimiter);
    	int [] uniqueGeneCounts= new int [genes.length];
    	for (int i=0; i<genes.length; i++) {
    		String g = genes[i];
    		Integer count = uniqueGeneCounter.getCountForKey(g);
    		if (count==null)
				count = new Integer (0);
    	    uniqueGeneCounts[i]=count;
    	}
    	return uniqueGeneCounts;
	}

	// for each primary read with 0 or more secondary reads, take all the gene tags and sort them in alphabetical order.
	// Make them into a single colon separated string and concatenate them, then increment the number of primary reads.
	/**
	 * For a list of reads with the same read name,
	 * @param records
	 * @param geneClusterCounter
	 * @return
	 */
	private ObjectCounter<String> populateMultiMapperGenes (final List<SAMRecord> records, final ObjectCounter<String> geneClusterCounter) {

		// iterate through the reads at the correct quality, looking for multi-mapper reads.
		Iterator<SAMRecord> filteredReads =  new MapQualityFilteringIterator(records.iterator(), this.MIN_READ_MQ, this.MAX_READ_MQ);
		filteredReads = new EditDistanceFilteringIterator(filteredReads, this.MAXIMUM_READ_EDIT_DISTANCE);
		List<String> geneList = getGeneList(filteredReads);

		// you need to see at least 2 genes to have an interesting multi-mapper.
		// mapping to a gene and a non-gene qualifies.
		if (geneList.size()<2) return (geneClusterCounter);
		String result = StringUtil.join(this.DELIMITER, geneList);
		geneClusterCounter.increment(result);
		return geneClusterCounter;
	}

	private ObjectCounter<String> populateUniqueMapperGenes (final List<SAMRecord> records, final ObjectCounter<String> uniqueGeneCounter) {
		Iterator<SAMRecord> filteredReads = new MapQualityFilteredIterator(records.iterator(), this.UNIQUE_READ_MQ, false);

		// Iterator<SAMRecord> filteredReads =  new MapQualityFilteringIterator(records.iterator(), this.UNIQUE_READ_MQ, this.UNIQUE_READ_MQ);

		List<String> geneList = getGeneList(filteredReads);
		if (geneList.size()>1)
			log.warn("More than 1 gene for a unique read?  Impossible?");

		for (String s: geneList)
			uniqueGeneCounter.increment(s);
		return uniqueGeneCounter;
	}

	/**
	 * For an iterator of reads, get an alphabetically ordered list of unique genes.
	 * If the iterator contains a read that does not overlap a gene, then the list will contain the
	 * MISSING_GENE string.
	 * TODO: Reads should be on the same strand as the annotation.  This can change via the strand strategy later
	 * if the newer version of the geneExonTagger is brought into the main codebase.
	 * @param iter
	 * @return
	 */
	private List<String> getGeneList (final Iterator<SAMRecord> iter) {
		Set<String> geneExonTagList = new HashSet<>();
		while (iter.hasNext()) {
			SAMRecord r = iter.next();
			String gene = r.getStringAttribute(this.GENE_EXON_TAG);
			String geneStrand = r.getStringAttribute(this.STRAND_TAG);
			if (gene==null) gene=MISSING_GENE;
			// if the strand isn't null, assert that the gene tag strand and read strand are the same.
			if (geneStrand!=null) {
				String readStrandString = Utils.strandToString(!r.getReadNegativeStrandFlag());
				if (geneStrand.equals(readStrandString))
					geneExonTagList.add(gene);
			} else // if the strand is null, the gene is too, so add it since it's the MISSING_GENE.
				geneExonTagList.add(gene);

		}
		List<String> geneList = new ArrayList<>(geneExonTagList);
		Collections.sort(geneList);
		return (geneList);
	}

	private GroupingIterator<SAMRecord> getGroupedReadNameIterator (final SamReader reader) {
		Iterator<SAMRecord> iter = CustomBAMIterators.getQuerynameSortedRecords(reader);
		GroupingIterator<SAMRecord> gIter = new GroupingIterator<>(iter, new ReadNameComparator());
		return (gIter);
	}

	/**
	 * A faster iterator that throws out more reads by filtering to a minimum read quality before queryname sorting.
	 * Queryname sorting is ignored if the BAM is already in queryname order.
	 * This iterator does not return all reads in the BAM.
	 * @param reader
	 * @return An iterator of reads with that have a minimum map quality of at least MIN_READ_MQ, in queryname order.
	 */
	@SuppressWarnings("resource")
	private GroupingIterator<SAMRecord> getGroupedReadNameIteratorFaster (final SamReader reader) {

		// MapQualityFilteredIterator
		Iterator<SAMRecord> iter = new MapQualityFilteredIterator(reader.iterator(), this.MIN_READ_MQ, false);

		// if the SamReader is NOT in queryName order, then wrap the iterator in a query name ordering iter.
		if (!reader.getFileHeader().getSortOrder().equals(SortOrder.queryname))
			iter = SamRecordSortingIteratorFactory.create(reader.getFileHeader(), iter, new SAMRecordQueryNameComparator(), this.pl);

		GroupingIterator<SAMRecord> gIter = new GroupingIterator<>(iter, new ReadNameComparator());
		return (gIter);
	}

	private class MapQualityFilteringIterator extends FilteredIterator<SAMRecord> {
		private int minMapQuality;
		private int maxMapQuality;
		public MapQualityFilteringIterator (final Iterator<SAMRecord> underlyingIterator, final int minMapQuality, final int maxMapQuality) {
			super(underlyingIterator);
			this.minMapQuality=minMapQuality;
			this.maxMapQuality=maxMapQuality;
		}

		@Override
	    protected boolean filterOut(final SAMRecord r) {
			int mq= r.getMappingQuality();
	    	if (mq<=this.maxMapQuality && mq>=this.minMapQuality) return false;
	    	return true;
	    }
	}

	private class EditDistanceFilteringIterator extends FilteredIterator<SAMRecord> {
		private Integer maxEditDistance;
		private static final String EDIT_DISTANCE_TAG="NM";

		public EditDistanceFilteringIterator (final Iterator<SAMRecord> underlyingIterator, final Integer maxEditDistance) {
			super(underlyingIterator);
			this.maxEditDistance=maxEditDistance;
		}

		@Override
	    protected boolean filterOut(final SAMRecord r) {
			if (maxEditDistance==null) return false;
			Object o = r.getAttribute(EDIT_DISTANCE_TAG);
			if (o==null) return false;
			int readEditDistance = ((Integer) o).intValue();
	    	return readEditDistance > this.maxEditDistance;
	    }


	}


	private List<MetaGene> buildMetaGenes (final ObjectCounter<String> geneClusterCounter, final ObjectCounter<String> uniqueGeneCounter) {
		StringInterner i = new StringInterner();
		List<MetaGene> result = new ArrayList<>();

		List<String> keys = geneClusterCounter.getKeysOrderedByCount(true);
        for (String k: keys) {
        	int countMetaGene = geneClusterCounter.getCountForKey(k);
        	String [] genes = StringUtils.split(k, this.DELIMITER);
        	int [] uniqueGeneCounts=getUniqueGeneCounts(k, this.DELIMITER, uniqueGeneCounter);
        	// intern all the strings!
        	for (int index=0; index<genes.length; index++)
				genes[index]=i.intern(genes[index]);
        	i.intern(k);
        	MetaGene metaGene = new MetaGene(i.intern(k), countMetaGene, genes, uniqueGeneCounts, this.MISSING_GENE);
        	result.add(metaGene);
        }

        Collections.sort(result, TagMetaGenes.PCT_MULTIMAP_NAME);
		return (result);

	}
	private class MetaGene {

		private String metaGeneName;
		private int metaGeneReadCount;
		private String [] uniqueGeneNames;
		private int [] uniqueGeneCounts;
		private int [] nonMissingGeneIndex;

		public MetaGene (final String metaGeneName, final int metaGeneReadCount, final String [] uniqueGenes, final int [] uniqueCounts, final String missingGeneString) {
			this.metaGeneName=metaGeneName;
			this.metaGeneReadCount=metaGeneReadCount;
			this.uniqueGeneNames=uniqueGenes;
			this.uniqueGeneCounts=uniqueCounts;

			// find missing gene positions.
			List<Integer> nonMissingGeneIndex = new ArrayList<>();
			for (int i=0; i<uniqueGeneCounts.length; i++)
				if (uniqueGeneNames[i].equals(missingGeneString)==false)
					nonMissingGeneIndex.add(i);

			// annoying that I can't autobox this.
			this.nonMissingGeneIndex = ArrayUtils.toPrimitive(nonMissingGeneIndex.toArray(new Integer[nonMissingGeneIndex.size()]));
		}

		public String [] getUniqueGeneNames (final boolean filterNonGenes) {
			if (filterNonGenes) {
				String [] result = new String [this.nonMissingGeneIndex.length];
				for (int i=0; i<this.nonMissingGeneIndex.length; i++)
					result[i]=this.uniqueGeneNames[this.nonMissingGeneIndex[i]];
				return (result);
			}
			return this.uniqueGeneNames;
		}

		public int [] getUniqueGeneCounts (final boolean filterNonGenes) {
			if (filterNonGenes) {
				int  [] result = new int [this.nonMissingGeneIndex.length];
				for (int i=0; i<this.nonMissingGeneIndex.length; i++)
					result[i]=this.uniqueGeneCounts[this.nonMissingGeneIndex[i]];
				return (result);
			}
			return this.uniqueGeneCounts;
		}

		public String getMetaGeneName() {
			return metaGeneName;
		}

		public int getMetaGeneReadCount() {
			return metaGeneReadCount;
		}

		public int getTotalUniqueReads (final boolean filterNonGenes) {
			int [] counts = getUniqueGeneCounts(filterNonGenes);
			int sumUGC = IntStream.of(counts).sum();
			return sumUGC;
		}

		public double getPercentMultiMappedReads (final boolean filterNonGenes) {
			int sumUGC=getTotalUniqueReads(filterNonGenes);
			double pctUnique = ((double) this.metaGeneReadCount / (sumUGC+this.metaGeneReadCount))*100;
			return pctUnique;
		}

		public boolean hasReadsOnNonGenicRegions () {
			return this.nonMissingGeneIndex.length!=this.uniqueGeneNames.length;
		}

	}

	/**
	 * Provides sorting of MetaGene objects first by their percent multimapping reads (excluding reads that don't align to a gene), then by name.
	 */
	private static final Comparator<MetaGene> PCT_MULTIMAP_NAME = new Comparator<MetaGene>() {
		@Override
		public int compare(final MetaGene o1, final MetaGene o2) {
			double p1=MathUtil.round(o1.getPercentMultiMappedReads(true), 1);
			double p2=MathUtil.round(o2.getPercentMultiMappedReads(true), 1);
		return new CompareToBuilder()
				.append(p2, p1)
				.append(o2.getMetaGeneReadCount(), o1.getMetaGeneReadCount())
				.append(o1.metaGeneName, o2.metaGeneName)
				.toComparison();
		}
	};


	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new TagMetaGenes().instanceMain(args));
	}
}


