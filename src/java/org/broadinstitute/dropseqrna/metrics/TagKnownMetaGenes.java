package org.broadinstitute.dropseqrna.metrics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.lang.StringUtils;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.SamHeaderUtil;

import picard.cmdline.CommandLineProgram;
import picard.cmdline.CommandLineProgramProperties;
import picard.cmdline.Option;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.BasicInputParser;

@CommandLineProgramProperties(
        usage = "A special case tagger.  Adds a new tag for meta-genes, which are genes known to have high homology, where relatively few reads will map uniquely to these regions.  Instead,"
        		+ "reads that overlap these regions receive a new tag, and all the reads under the set of genes represented by this tag's reads can be pooled under other analysis like DGE. "
        		+ "These meta-genes are determined by the user, who passes in a specified list of genes that should be aggregated." ,
        usageShort = "Tags meta-genes",
        programGroup = DropSeq.class
)
public class TagKnownMetaGenes extends CommandLineProgram {

	private final Log log = Log.getInstance(TagReadWithGeneExon.class);
	private ProgressLogger pl = new ProgressLogger(log);

	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze")
	public File INPUT;

	@Option(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output BAM, written with new Gene/Exon tag")
	public File OUTPUT;

	@Option(doc="A file containing a list of meta gene names associated to one or more genes.  The first column contains a metagene name, "
			+ "and the other column contains a comma separated list of gene names")
	public File METAGENES;

	@Option(doc = "The tag to combine.")
	public String TAG="GE";

	@Option (doc="The metagene tag to add to records")
	public String METAGENE_TAG="MG";

	@Option(doc="The minimum map quality of the read to get a metatag.")
	public Integer MIN_READ_MQ=3;

	@Option(doc="The maximum map quality of the read to get a metatag.")
	public Integer MAX_READ_MQ=3;

	@Override
	protected int doWork() {
		IOUtil.assertFileIsReadable(this.INPUT);
		IOUtil.assertFileIsReadable(this.METAGENES);
		IOUtil.assertFileIsWritable(this.OUTPUT);
		Map<String, String> metaGeneMap = readMetaGeneFile(METAGENES);

		SamReader inputSam = SamReaderFactory.makeDefault().open(INPUT);

		SAMFileHeader header = inputSam.getFileHeader();
		SamHeaderUtil.addPgRecord(header, this);
		SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT);

		for (SAMRecord r: inputSam) {
        	pl.record(r);

        	if (testMapQuality(r, this.MIN_READ_MQ, this.MAX_READ_MQ))
				r=	setMetaGenes(r, metaGeneMap);
        	writer.addAlignment(r);
        }

		CloserUtil.close(inputSam);
		writer.close();

		return 0;
	}

	/**
	 * Returns false for reads that aren't between the min and max map quality thresholds, or are unmapped.
	 * @param r
	 * @param minMapQual
	 * @param maxMapQual
	 * @return
	 */
	private boolean testMapQuality (final SAMRecord r, final int minMapQual, final int maxMapQual) {
		if (r.getReadUnmappedFlag()) return false;
		int mq = r.getMappingQuality();
		if (mq>=minMapQual && mq <=maxMapQual) return true;
		return false;
	}

	private SAMRecord setMetaGenes (final SAMRecord r, final Map<String, String> metaGeneMap) {
		Object tagValue = r.getAttribute(this.TAG);
		if (tagValue==null) return (r);
		String tv = (String) tagValue;

		// look up gene in list.
		String metaGene = metaGeneMap.get(tv);
		if (metaGene==null) return (r); // the read doesn't contain a tag that needs a meta tag.
		// read has a recognized tag, add it.
		r.setAttribute(this.METAGENE_TAG, metaGene);
		return r;
	}


	/**
	 * Parse a meta-gene file.  First column has the meta-gene name.
	 * Second column is a comma separated list of genes.
	 * This returns a map with a key of the gene name, and a value of the meta-gene name.
	 * @param input
	 * @return
	 */
	public static Map<String, String> readMetaGeneFile (final File input) {
		IOUtil.assertFileIsReadable(input);
		Map<String, String> result = new HashMap<String,String>();
		BasicInputParser parser = new BasicInputParser(false, 2, input);
		// validate header.
		if (parser.hasNext()) {
			String [] line =parser.next();
			if (!line[0].equals("METAGENE") || !line[1].equals("GENES")) {
				parser.close();
				throw new IllegalArgumentException("Wrong header, is this the metagene file you're looking for? " + parser.getCurrentLine());
			}
		}
		while(parser.hasNext()) {
			String [] line =parser.next();
			String [] genes = StringUtils.split(line[1], ',');
			for (String g: genes)
				result.put(g, line[0]);
		}
		parser.close();
		return (result);
	}

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new TagKnownMetaGenes().instanceMain(args));
	}
}
