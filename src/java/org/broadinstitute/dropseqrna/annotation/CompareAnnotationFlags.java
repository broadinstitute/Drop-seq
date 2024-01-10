package org.broadinstitute.dropseqrna.annotation;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.annotation.functionaldata.*;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.readiterators.GeneFunctionProcessor;
import org.broadinstitute.dropseqrna.utils.readiterators.StrandStrategy;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.samtools.util.StringUtil;
import picard.annotation.LocusFunction;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

@CommandLineProgramProperties(
        summary = "Test the old vs new versions of TagReadWithGeneFunction by comparing the tags from both for UTR/CODING reads.",
        		oneLineSummary ="Test program don't use.",
        programGroup = DropSeq.class)

@groovy.transform.Generated
public class CompareAnnotationFlags extends CommandLineProgram {


	private final Log log = Log.getInstance(CompareAnnotationFlags.class);
	private ProgressLogger pl = new ProgressLogger(log);


	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input SAM or BAM file to analyze")
	public File INPUT;

	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc = "The output BAM, has reads that are different between the two taggers.")
	public File OUTPUT;

	@Argument(doc="Output BAM for split-read bugs.")
	public File OUTPUT_SPLIT_READ_BUG;

	@Argument (doc="Combinations of genes that lead to ambiguous assignment.")
	public File AMBIGUOUS_GENE_OUTPUT;

	@Argument(doc="Gene Exon tag.  Old version of TagReadWithGeneExon tag.")
	public String GE_TAG="GE";

	@Argument(doc="Gene Stand tag.  Old version of TagReadWithGeneExon tag")
	public String GS_TAG="GS";

	@Argument(doc="Gene Function tag.  Old version of TagReadWithGeneExon tag")
	public String FUNCTION_TAG="XF";

	@Argument(doc="Gene Name tag.  Takes on the gene name this read overlaps (if any)")
	public String GENE_NAME_TAG="gn";

	@Argument(doc="Gene Strand tag.  For a given gene name <GENE_NAME_TAG>, this is the strand of the gene.")
	public String GENE_STRAND_TAG="gs";

	@Argument(doc="Gene Function tag.  For a given gene name <GENE_NAME_TAG>, this is the function of the gene at this read's position: UTR/CODING/INTRONIC/...")
	public String GENE_FUNCTION_TAG="gf";

	private String DELIMITER = ",";

	private String [] problemGenes={"CDK11A", "RP1-283E3.8"};

	private int splitReadBugCount=0;
	private int problemGeneReadCount=0;

	// @Argument (doc="A list of functional annotations that reads need to be completely contained by to be considered for analysis.")
    public List<LocusFunction> LOCUS_FUNCTION_LIST=new ArrayList<LocusFunction>(Arrays.asList(LocusFunction.CODING, LocusFunction.UTR));

	// genes that might have the wrong results because of the split-read bug.
	private Set<String> newGenesHigher = new HashSet<String>();
	private Set<String> oldGenesHigher = new HashSet<String>();
	private Set<String> totalGenes = new HashSet<String>();
	private int counterOldReadHigher=0;
	private int counterNewReadHigher=0;
	private int totalReads =0;

	@Override
	protected int doWork() {
		IOUtil.assertFileIsReadable(this.INPUT);
		IOUtil.assertFileIsWritable(this.OUTPUT);
		IOUtil.assertFileIsWritable(this.OUTPUT_SPLIT_READ_BUG);
		if (AMBIGUOUS_GENE_OUTPUT!=null) IOUtil.assertFileIsWritable(AMBIGUOUS_GENE_OUTPUT);

		SamReader inputSam = SamReaderFactory.makeDefault().open(INPUT);
		SAMFileHeader header = inputSam.getFileHeader();
		SAMFileWriter writer= new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT);
		SAMFileWriter writer2= new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT_SPLIT_READ_BUG);

		int counter=0;
		int totalNumReads=0;
		LocusFunction [] acceptedLocusFunctions = {LocusFunction.UTR, LocusFunction.CODING};
		GeneFunctionProcessor gfiw = new GeneFunctionProcessor("gn", "gs", "gf", false, StrandStrategy.SENSE, Arrays.asList(acceptedLocusFunctions));

		FunctionalDataProcessorI fdp =  FunctionalDataProcessorFactory.getFunctionalDataProcessor(StrandStrategy.SENSE, LOCUS_FUNCTION_LIST, FunctionalDataProcessorStrategyEnum.DROPSEQ);
		ObjectCounter<String> ambiguousGeneCounter = new ObjectCounter<String>();

		for (SAMRecord r: inputSam) {
			pl.record(r);
			totalNumReads++;
			boolean flag = testTagsDisagree(r, gfiw, writer2);

			// test and record ambiguous genes.
			Set<String> ambiguousGenes = findAmbiguousGenes(r, fdp);
			if (ambiguousGenes.size()>0)
				ambiguousGeneCounter.increment(orderAndContatonate(ambiguousGenes));
			if (flag) {
				writer.addAlignment(r);
				counter++;
			}


		}
		if (this.AMBIGUOUS_GENE_OUTPUT!=null) writeAmbiguousGenes(this.AMBIGUOUS_GENE_OUTPUT, ambiguousGeneCounter);

		writer.close();
		writer2.close();
		log.info("Number of reads: " + totalNumReads +" Number with disagreements " + counter + " Split read bug reads " + this.splitReadBugCount + " Problem Gene reads " + this.problemGeneReadCount);
		log.info("#Genes with higher expression old data:" + this.oldGenesHigher.size());
		log.info("#Genes with higher expression new data:" + this.newGenesHigher.size());
		log.info("#Total Genes " + this.totalGenes.size());
		log.info("Number of old reads rejected: " + this.counterOldReadHigher);
		log.info("Number of new reads rejected: " + this.counterNewReadHigher);
		log.info("Total # reads with expression: " + this.totalReads);
		return 0;
	}

	private String orderAndContatonate (final Collection<String> strings) {
		StringBuilder result = new StringBuilder();
		List<String> l = new ArrayList<String>(strings);
		Collections.sort(l);

		Iterator<String> iter = l.iterator();
		if (!iter.hasNext()) return result.toString();
		result.append(iter.next());

		while (iter.hasNext()) {
			result.append(",");
			result.append(iter.next());
		}
		return result.toString();
	}
	/**
	 * For reach read, check if it has coding on multiple genes.
	 * Build up a count of what genes are involved.  Are certain genes (or patterns of genes like CTD or RP) common?
     * 	1) Order the annotation types to pick the "preferred" one (UTR/coding/intronic).
	 *  2) Find reads with > 1 gene at the preferred annotation type.
	 * @param r
	 * @return
	 */
	private Set<String> findAmbiguousGenes (final SAMRecord r , final FunctionalDataProcessorI fdp) {
		Set<String> result = new HashSet<String>();


		String geneList = r.getStringAttribute(this.GENE_NAME_TAG);
		String strandList = r.getStringAttribute(this.GENE_STRAND_TAG);
		String functionList = r.getStringAttribute(this.GENE_FUNCTION_TAG);
		// If you're missing the gene, strand, or function, you can't use this
				// read.
		if (geneList == null || strandList == null || functionList == null)
			return result;
		String[] genes = geneList.split(DELIMITER);
		String[] strands = r.getStringAttribute(GENE_STRAND_TAG).split(DELIMITER);
		LocusFunction[] locusFunctions = getLocusFunctionFromRead(functionList);

		List<FunctionalData> fdList = fdp.getFilteredFunctionalData(genes, strands, locusFunctions, r.getReadNegativeStrandFlag());
		fdList = fdp.filterToPreferredAnnotations(fdList);

		// if you have > 1 assignment, there's ambiguity.
		if (fdList.size()>1)
			for (FunctionalData fd: fdList)
				result.add(fd.getGene());
		return result;
	}

	private static void writeAmbiguousGenes(final File file, final ObjectCounter<String> ambiguousGeneCounter) {
        final BufferedWriter writer = IOUtil.openFileForBufferedWriting(file);
        List<String> genes = ambiguousGeneCounter.getKeysOrderedByCount(true);
        try {
        	String [] header = {"GENES", "COUNT"};
        	writer.write(StringUtil.join("\t", header));
        	writer.newLine();
            for (final String gene : genes) {
            	String [] line = {gene, Integer.toString(ambiguousGeneCounter.getCountForKey(gene))};
            	String result = StringUtil.join("\t", line);
                writer.write(result);
                writer.newLine();
            }
            writer.close();
        } catch (IOException e) {
            throw new RuntimeIOException("Exception writing " + file.getAbsolutePath(), e);
        }
    }

	/**
	 * If
	 * @param r
	 * @return
	 */
	private boolean testTagsDisagree (final SAMRecord r, final GeneFunctionProcessor gfiw, final SAMFileWriter splitReadWriter) {
		// for the old data, if there's expression oldGenes size =1.
		Set<String> oldGenes = getCodingGenesOld(r);
		Set<String> newGenes = getCodingGenesNew(r);
		totalGenes.addAll(oldGenes);
		totalGenes.addAll(newGenes);

		// For the new data, if there's expression, #recs =1.  Otherwise recs=0.
		List<SAMRecord> recs = gfiw.processRead(r);

		boolean oldExpression = oldGenes.size()==1;
		boolean newExpression = recs.size()==1;

		if (oldExpression!=newExpression) {
			if (oldExpression && oldGenes.size()> recs.size()) {
				this.oldGenesHigher.addAll(oldGenes);
				this.counterOldReadHigher++;
			}
			if (newExpression && recs.size()> oldGenes.size()) {
				String gene = recs.get(0).getStringAttribute("gn");
				this.newGenesHigher.add(gene);
				this.counterNewReadHigher++;
			}


		}
		if (oldGenes.size()==0 & newGenes.size()==0) return false;
		totalReads++;

		// if the new tagger contains multiple coding genes on the same strand as the read, then we expect the old tagger to have no genes tagged.
		if (newGenes.size()>1 & oldGenes.size()==0) return false;

		// if the genes are in the problem genes set, ignore this read's problem.
		if (testProblemGenes(oldGenes, newGenes)) return false;

		if (testSplitReadBug(r)) {
			splitReadWriter.addAlignment(r);
			return false;
		}

		// do we have the same set of genes?
		boolean flag = oldGenes.containsAll(newGenes) && newGenes.containsAll(oldGenes);
		/*
		if (!flag)
			log.info("Rejected read" + r.getReadName());
		*/

		return !flag;

	}

	/**
	 * Is the gene tagged on the whitelist of problem genes?
	 * @param oldGenes
	 * @param newGenes
	 * @return
	 */
	private boolean testProblemGenes (final Set<String> oldGenes, final Set<String> newGenes) {
		for (String g: this.problemGenes)
			if (oldGenes.contains(g) || newGenes.contains(g)) {
				problemGeneReadCount++;
				return true;
			}
		return false;
	}


	/**
	 * Find cases where the old tagger gave no tag, but the new tagger provides both an exonic and an intronic tag on the same gene, and the read is split across exon/exon(?) boundaries.
	 * This should find any gene with both an intron and coding tag for the same gene, and remove those genes from the new tagger.  Then compare the new and old gene sets remaining.
	 * @param r
	 * @return true if this read is different because of the split-read bug.
	 */
	private boolean testSplitReadBug(final SAMRecord r) {
		if (r.getReadName().equals("000000000-AMY9M:1:1104:11511:11301"))
			log.info("STOP");
		Set<String> oldGenes = getCodingGenesOld(r);
		Set<String> newGenes = getCodingGenesNew(r);
		// gene sets match up already.
		if (oldGenes.containsAll(newGenes) && newGenes.containsAll(oldGenes))
			return false;

		Set<String> buggedGenes = splitReadBugGenes(r);
		ObjectCounter<Integer> conditions = new ObjectCounter<Integer>();
		// if a single gene appears in the old set, and that gene appears in the bugged gene list and the newGenes has > 1 gene, then there will be a change in expression (old will be higher)
		if (oldGenes.size()==1 & newGenes.size()>1)
			conditions.increment(1);
			// this.oldGenesHigher.addAll(oldGenes);
		// if no gene appears in the old set and greater than one gene appears in the new set then expression shouldn't change
		if (oldGenes.size()==0 & newGenes.size()>1)
			conditions.increment(2);

		// if no gene appears in the old set and one gene appears in the new set, then expression will change (new will be higher)
		if (oldGenes.size()==0 & newGenes.size()==1 & buggedGenes.containsAll(newGenes))
			conditions.increment(3);
			// newGenesHigher.addAll(newGenes);
		if (conditions.getSize()>1)
			log.info(r.getReadName()+  " Multiple conditions");


		newGenes.removeAll(buggedGenes);


		if (oldGenes.containsAll(newGenes) & buggedGenes.size()>0) {
			this.splitReadBugCount++;
			return true;
		}

		return false;
	}

	private Set<String> splitReadBugGenes (final SAMRecord r) {
		Set<String> result = new HashSet<String>();

		String geneList = r.getStringAttribute(this.GENE_NAME_TAG);
		String strandList = r.getStringAttribute(this.GENE_STRAND_TAG);
		String functionList = r.getStringAttribute(this.GENE_FUNCTION_TAG);
		if (geneList==null || strandList==null || functionList==null) return result;

		String[] genes = geneList.split(DELIMITER);
		String[] strands = strandList.split(DELIMITER);
		LocusFunction[] locusFunctions = getLocusFunctionFromRead(functionList);

		Set<String> geneSet = new HashSet<String>();
		geneSet.addAll(Arrays.asList(genes));
		for (String g: geneSet)
			if (hasSplitReadBug(g, genes, strands, locusFunctions, r))
				result.add(g);
		return result;
	}

	private boolean hasSplitReadBug (final String gene, final String[] genes, final String[] strands, final LocusFunction[] locusFunctions, final SAMRecord r) {
		// if the read isn't split, this bug will not occur.
		if (!readIsSplit(r)) return false;

		ObjectCounter<LocusFunction> locusFuncCounter = new ObjectCounter<LocusFunction>();
		for (int i=0; i<genes.length; i++) {
			boolean correctStrand=(strands[i].equals("+") && !r.getReadNegativeStrandFlag()) || (strands[i].equals("-") && r.getReadNegativeStrandFlag());
			if (genes[i].equals(gene) && correctStrand)
				locusFuncCounter.increment(locusFunctions[i]);
		}

		// should have two locus functions.
		int countCoding = locusFuncCounter.getCountForKey(LocusFunction.CODING);
		int countUTR = locusFuncCounter.getCountForKey(LocusFunction.UTR);
		int countIntronic = locusFuncCounter.getCountForKey(LocusFunction.INTRONIC);

		// read should be split by cigar.



		if ((countCoding==1 || countUTR==1) && countIntronic==1) {
			this.splitReadBugCount++;
			return true;
		}
		return false;
	}

	private boolean readIsSplit (final SAMRecord r) {
		List<CigarElement> cigarElements = r.getCigar().getCigarElements();
		int countM =0;
		for (CigarElement e: cigarElements)
			if (e.getOperator()==CigarOperator.M)
				countM++;
		if (countM>1) return true;
		return false;
	}

	private Set<String> getCodingGenesOld (final SAMRecord r) {
		Set<String> result = new HashSet<String>();

		String oldGene = r.getStringAttribute(this.GE_TAG);
		String oldFunction = r.getStringAttribute(this.FUNCTION_TAG);
		String oldStrand = r.getStringAttribute(this.GS_TAG);
		if (oldGene==null || oldStrand==null || oldFunction==null) return result;
		LocusFunction lf = LocusFunction.valueOf(oldFunction);

		boolean hasCoding=lf==LocusFunction.CODING || lf==LocusFunction.UTR;
		// is the read on the same strand?
		boolean correctStrand = (oldStrand.equals("+") && !r.getReadNegativeStrandFlag()) || (oldStrand.equals("-") && r.getReadNegativeStrandFlag());
		if (hasCoding && correctStrand) result.add(oldGene);
		return (result);

	}

	private Set<String> getCodingGenesNew (final SAMRecord r) {
		Set<String> result = new HashSet<String>();

		String geneList = r.getStringAttribute(this.GENE_NAME_TAG);
		String strandList = r.getStringAttribute(this.GENE_STRAND_TAG);
		String functionList = r.getStringAttribute(this.GENE_FUNCTION_TAG);

		// if anything is null, return empty list.
		if (geneList==null || strandList==null || functionList==null) return result;

		String[] genes = geneList.split(DELIMITER);
		String[] strands = strandList.split(DELIMITER);
		LocusFunction[] locusFunctions = getLocusFunctionFromRead(functionList);

		for (int i=0; i<genes.length; i++) {
			boolean hasCoding= (locusFunctions[i]==LocusFunction.CODING || locusFunctions[i]==LocusFunction.UTR);
			boolean correctStrand=(strands[i].equals("+") && !r.getReadNegativeStrandFlag()) || (strands[i].equals("-") && r.getReadNegativeStrandFlag());
			if (hasCoding && correctStrand)
				result.add(genes[i]);
		}
		return (result);

	}

	private LocusFunction[] getLocusFunctionFromRead(final String functionList) {
		String[] fl = functionList.split(DELIMITER);
		LocusFunction[] result = new LocusFunction[fl.length];
		for (int i = 0; i < fl.length; i++) {
			LocusFunction lf = LocusFunction.valueOf(fl[i]);
			result[i] = lf;
		}
		return result;
	}

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new CompareAnnotationFlags().instanceMain(args));
	}

}
