/*
 * The Broad Institute
 * SOFTWARE COPYRIGHT NOTICE AGREEMENT
 * This software and its documentation are copyright 2017 by the
 * Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
 *
 * This software is supplied without any warranty or guaranteed support whatsoever.
 * Neither the Broad Institute nor MIT can be responsible for its use, misuse,
 * or functionality.
 */

package org.broadinstitute.dropseqrna.eqtl;

import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.*;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.rank.Percentile;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.annotation.GeneAnnotationReader;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.tools.DGEMatrix;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.tools.MatrixTransformFactory;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.tools.MatrixTransformI;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.IntervalTagComparator;
import org.la4j.Matrix;
import org.la4j.Vectors;
import org.la4j.vector.functor.VectorFunction;
import picard.annotation.Gene;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        summary = "Generate expression matrix data for Matrix EQTL and tensorQTL packages",
        oneLineSummary = "Generate expression matrix data for Matrix EQTL and tensorQTL packages",
        programGroup = DropSeq.class
)
public class PrepareEqtlExpressionData extends CommandLineProgram {

	private final Log log = Log.getInstance(PrepareEqtlExpressionData.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input meta cell expression files.")
	public List<File> META_CELL_FILE;

	@Argument (doc="A file with 1 column no header, containing identifiers of donors that should be excluded from the expression output.", optional=true)
	public File REJECTED_DONOR_LIST;

	@Argument (doc = "The annotations file that provides gene locations.  Supports GTF and RefFlat format.")
	public File ANNOTATIONS_FILE;

	@Argument (doc="Sequence dictionary file")
	public File SD_FILE;

	@Argument(
			doc = "The output expression matrix.  First entry is the gene symbol, columns following that "
					+ "are each donor's expression.  One gene per row.",
			optional = true
	)
	public File EXPRESSION_FILE;

	@Argument (doc="Locations of the genes in the expression matrix", optional = true)
	public File GENE_LOCATION_FILE;

	@Argument(
			doc="The output phenotype expression BED. "
					+ "See tensorQTL for more details: " 
					+ "https://github.com/broadinstitute/tensorqtl/tree/master#input-formats",
			optional = true
	)
	public File EXPRESSION_BED_FILE;

	@Argument (doc="A list of donors which may contain repeats.  All donors from the META_CELL_FILES, with the rejected donors from the REJECTED_DONOR_LIST removed.")
	public File OUT_DONOR_LIST;

	@Argument (doc="Each meta cell will have this many total transcripts distributed across all genes.  Set this value to -1 to not normalize donor expression to be"
			+ "new expression = expression gene/sum(all genes) *TRANSCRIPTS_PER_CELL ")
	public Integer TRANSCRIPTS_PER_CELL=100000;

	@Argument (doc="Remove the bottom PCT of expressed genes from the data set. If unset no data is removed.  This calculates the median expression across for each gene "
			+ "across meta cells, then orders genes by expression and calculates percentiles.  The bottom <REMOVE_PCT_EXPRESSION> expressed genes are filtered.", optional=true)
	public Double REMOVE_PCT_EXPRESSION=50d;

	@Argument (doc="If true, donors with the same names from different experiments have their UMIs combined.")
	public boolean MERGE_DONORS=false;

	private final String FIELD_SEPARATOR="\t";

	@Argument (doc="A list of chromosomes to omit from the analysis.  The default is to omit the Y chromosome and MT.")
	public List<String> IGNORED_CHROMOSOMES= new ArrayList<>(Arrays.asList("Y", "MT"));

	private static final NumberFormat EXPRESSION_FORMATTER = new DecimalFormat("0.#####");

	@Override
	public int doWork() {
		if (META_CELL_FILE.size()==0) {
			log.info("Must supply at least one expression meta cell file");
			return 1;
		}
		validateMetaCellFiles(this.META_CELL_FILE);
		IOUtil.assertFileIsReadable(ANNOTATIONS_FILE);
		IOUtil.assertFileIsReadable(SD_FILE);

		if (EXPRESSION_FILE != null && GENE_LOCATION_FILE == null) {
			log.info(
					"If EXPRESSION_FILE is supplied, GENE_LOCATION_FILE must be supplied."
			);
			return 1;
		}
		if (EXPRESSION_FILE == null && GENE_LOCATION_FILE != null) {
			log.info(
					"If GENE_LOCATION_FILE is supplied, EXPRESSION_FILE must be supplied."
			);
			return 1;
		}

		if (EXPRESSION_FILE != null) {
			IOUtil.assertFileIsWritable(EXPRESSION_FILE);
			IOUtil.assertFileIsWritable(GENE_LOCATION_FILE);
		}
		if (EXPRESSION_BED_FILE!=null) {
			IOUtil.assertFileIsWritable(EXPRESSION_BED_FILE);
		}
		IOUtil.assertFileIsWritable(OUT_DONOR_LIST);
		if (REJECTED_DONOR_LIST!=null) IOUtil.assertFileIsReadable(REJECTED_DONOR_LIST);

		processData(
				this.META_CELL_FILE, this.REJECTED_DONOR_LIST, this.ANNOTATIONS_FILE, this.SD_FILE,
				this.EXPRESSION_FILE, this.GENE_LOCATION_FILE, this.EXPRESSION_BED_FILE, this.OUT_DONOR_LIST,
				this.TRANSCRIPTS_PER_CELL, this.IGNORED_CHROMOSOMES, this.REMOVE_PCT_EXPRESSION
		);

		return 0;
	}

	/**
	 * The end to end processing of this class encapsulated in a single method to allow pipelining.
	 * See program options for details.
	 * @return A list of donor names, the same as the context of outDonorFile.
	 */
	public List<String> processData(
			final List<File> metaCellFiles, final File rejectedDonorListFile, final File annotationFile,
			final File sequenceDictionaryFile, final File expressionFile, final File geneLocationFile,
			final File expressionBedFile, final File outDonorFile, final int transcriptsPerCell,
			final List<String> ignoredChromosomes, final Double removePctExpression
	) {

		if (metaCellFiles.size()==0)
			throw new IllegalArgumentException("Must supply at least one expression meta cell file");
		validateMetaCellFiles(metaCellFiles);
		IOUtil.assertFileIsReadable(annotationFile);
		IOUtil.assertFileIsReadable(sequenceDictionaryFile);

		if (expressionFile != null && geneLocationFile == null) {
			throw new IllegalArgumentException(
					"If expressionFile is supplied, geneLocationFile must be supplied."
			);
		}
		if (expressionFile == null && geneLocationFile != null) {
			throw new IllegalArgumentException(
					"If geneLocationFile is supplied, expressionFile must be supplied."
			);
		}
		if (expressionFile == null && expressionBedFile == null) {
			throw new IllegalArgumentException(
					"At least one of expressionFile and expressionBedFile should be specified"
			);
		}
		if (expressionFile != null) {
			IOUtil.assertFileIsWritable(expressionFile);
			IOUtil.assertFileIsWritable(geneLocationFile);
		}
		if (expressionBedFile != null) {
			IOUtil.assertFileIsWritable(expressionBedFile);
		}

		IOUtil.assertFileIsWritable(outDonorFile);
		if (rejectedDonorListFile!=null) IOUtil.assertFileIsReadable(rejectedDonorListFile);

		Set<String> rejectedDonors = Collections.emptySet();

		if (rejectedDonorListFile!=null)
			rejectedDonors = new HashSet<> (ParseBarcodeFile.readCellBarcodeFile(rejectedDonorListFile));

		DGEMatrix dge = readAndFilterMatrix(metaCellFiles.get(0), rejectedDonors, 0+":");

		for (int i=1; i<metaCellFiles.size(); i++) {
			DGEMatrix other = readAndFilterMatrix(metaCellFiles.get(i), rejectedDonors, i+":");
			dge= dge.merge(other);
		}

		//filter genes to remove rejected contigs.
		OverlapDetector<Gene> geneOD = getGeneMap(annotationFile, sequenceDictionaryFile);
		geneOD=filterGeneOD(geneOD, ignoredChromosomes);

		//filter genes to the set that exist in the GTF file before writing expression matrix.
		dge = filterMatrixByGTF(dge, geneOD);

		// normalize expression of remaining genes.
		if (transcriptsPerCell!=-1) {
			dge.applyTransform(MatrixTransformFactory.normalizeColumns());
			dge.applyTransform(new Multiply(transcriptsPerCell));
		}
		
		// filter expression data to genes above median threshold.
		dge=filterByPctExpression(dge, removePctExpression);

		// get a list of genes that are in both the GTF and DGE, and in genomic order.
		List<Gene> orderedGeneList = getOrderedGeneList(dge, geneOD, sequenceDictionaryFile);

		if (expressionFile != null) {
			//write out the expression data.
			writeExpressionMatrix(dge, orderedGeneList, expressionFile);
			// write out the gene locations.
			writeGeneLocations(orderedGeneList, geneLocationFile);
		}

		if (expressionBedFile != null) {
			// write out the gene counts
			writeExpressionBed(dge, orderedGeneList, expressionBedFile);
		}

		// writeDonorNames, which are cleaned up to remove prefixes.
		List<String> donorNames = cleanupDonorNames(dge.getCellBarcodes());
		writeDonorNames(donorNames, outDonorFile);
		return donorNames;
	}

	private DGEMatrix filterByPctExpression (final DGEMatrix dge, final Double pct) {
		// no op.
		if (pct==null || pct==0) return dge;
		Mean m = new Mean();
		Map<String, Double> meanExpression = new HashMap<>();
		for (String gene: dge.getGenes()) {
			double [] exp = dge.getExpression(gene);
			double meanValue = m.evaluate(exp);
			meanExpression.put(gene, meanValue);
		}

		Percentile percentile = new Percentile();
		double [] expMeans = meanExpression.values().stream().mapToDouble(x -> x).toArray();
		double threshold = percentile.evaluate(expMeans, pct);

		List<String> genesToFilter = new ArrayList<>();
		for (String gene: meanExpression.keySet())
			if (meanExpression.get(gene) < threshold)
				genesToFilter.add(gene);

		dge.toDenseMatrix();
		dge.removeGenes(genesToFilter);

		return dge;
	}
	private OverlapDetector<Gene> filterGeneOD (final OverlapDetector<Gene> geneOD, final Collection <String> ignoredChromosomes) {
		Set<String> contigs = new HashSet<>(ignoredChromosomes);
		OverlapDetector<Gene> result = new OverlapDetector<>(0, 0);
		for (Gene g: geneOD.getAll())
			if (!contigs.contains(g.getContig()))
				result.addLhs(g, g);
		return result;
	}

	/**
	 * Need to build a final list of gene objects in genomic order, containing only the genes that are in the DGE.
	 * @param dge
	 * @param geneOD
	 * @return
	 */
	private List<Gene> getOrderedGeneList (final DGEMatrix dge, final OverlapDetector<Gene> geneOD, final File sequenceDictionaryFile) {
		// get the set of strings for searching.
		Set<String> dgeGenes = new HashSet<>(dge.getGenes());

		List<Gene> orderedGeneList = new ArrayList<>();
		for (Gene g: geneOD.getAll())
			if (dgeGenes.contains(g.getName()))
				orderedGeneList.add(g);
		// order them in genomic order, except that it treats contigs as strings.  Should use the sequence dictionary to properly sort contigs.
		// can I use IntervalTagComparator.compare to compare the intervals of the genes?
		// Collections.sort(orderedGeneList);
		SAMSequenceDictionary dict = getSD(sequenceDictionaryFile);
		GeneComparator gc = new GeneComparator(dict);
		orderedGeneList.sort(gc);
		return orderedGeneList;
	}

	private DGEMatrix filterMatrixByGTF (final DGEMatrix dge, final OverlapDetector<Gene> od) {
		dge.toDenseMatrix();
		// make a copy.
		Set<String> annoGenes = od.getAll().stream().map(Interval::getName).collect(Collectors.toSet());
		List<String> dgeGenes = new ArrayList<> (dge.getGenes());
		// dgeGenes has all genes that aren't in the annotations.
		dgeGenes.removeAll(annoGenes);
		dge.removeGenes(dgeGenes);
		return dge;
	}

	/**
	 * Load up the sequence dictionary and annotations file, and map gene names to gene objects.
	 * @param annotationsFile
	 * @param sequenceDictionary
	 * @return
	 */
	private OverlapDetector<Gene> getGeneMap (final File annotationsFile, final File sequenceDictionary) {
		SAMSequenceDictionary SD = getSD(sequenceDictionary);
		return GeneAnnotationReader.loadAnnotationsFile(annotationsFile, SD, VALIDATION_STRINGENCY);
	}

	private SAMSequenceDictionary getSD (final File sequenceDictionary) {
		SamReader in = SamReaderFactory.makeDefault().enable(SamReaderFactory.Option.EAGERLY_DECODE).open(sequenceDictionary);
		SAMSequenceDictionary result = in.getFileHeader().getSequenceDictionary();
		CloserUtil.close(in);
		return result;
	}

	/**
	 * Simple class to multiply the DGE values by a flat number.
	 * @author nemesh
	 *
	 */
	private class Multiply implements MatrixTransformI {
		private final double factor;
		public Multiply (final double factor) {
			this.factor=factor;
		}
		@Override
		public void apply(final Matrix m) {
			VectorFunction f = Vectors.asMulFunction(this.factor);
			int numRows = m.rows();
			 for (int i=0; i<numRows; i++)
				m.updateRow(i, f);
		}

	}

	private List<String> cleanupDonorNames (final Collection<String> donorNames) {
		return donorNames.stream().map(x -> x.split(":")[1]).collect(Collectors.toList());
	}

	private void writeDonorNames (final List<String> donorNames, final File output) {
		try {
			IOUtil.assertFileIsWritable(output);
			BufferedWriter out = new BufferedWriter(new OutputStreamWriter(IOUtil.openFileForWriting(output)));
			//write body.
			for (String name: donorNames) {
				out.write(name);
				out.newLine();
			}
			out.close();
		} catch (IOException e) {
			throw new RuntimeException("Trouble writing " + output.getAbsolutePath(), e);
		}
	}

	private void writeGeneLocations (final List<Gene> orderedGeneList, final File output) {
		try {
			IOUtil.assertFileIsWritable(output);
			BufferedWriter out = new BufferedWriter(new OutputStreamWriter(IOUtil.openFileForWriting(output)));
			//write header.
			String [] header = {"geneid", "chr", "s1", "s2"};
			out.write(StringUtils.join(header, this.FIELD_SEPARATOR));
			out.newLine();
			//write body.
			for (Gene g: orderedGeneList) {
				String [] body = {g.getName(), g.getContig(), Integer.toString(g.getStart()), Integer.toString(g.getEnd())};
				out.write(StringUtils.join(body, this.FIELD_SEPARATOR));
				out.newLine();
			}
			out.close();
		} catch (IOException e) {
			throw new RuntimeException("Trouble writing " + output.getAbsolutePath(), e);
		}
	}
	/**
	 * This is a hell of a lot like the DGE's writer, but has "id" instead of "gene" as the first row name.  It also skips the DgeHeader work that isn't needed.
	 * @param dge
	 * @param output
	 */
	private void writeExpressionMatrix(final DGEMatrix dge, final List<Gene> orderedGeneList, final File output) {

        try {
			IOUtil.assertFileIsWritable(output);
			List<String> cellBarcodes = dge.getCellBarcodes();

			BufferedWriter out = new BufferedWriter(new OutputStreamWriter(IOUtil.openFileForWriting(output)));

			// write header
			List<String> header = new ArrayList<>(cellBarcodes.size()+1);
			header.add("id");
			for (String c: cellBarcodes) {
				// get rid of the batch number :
				String cc = c.split(":")[1];
				header.add(cc);
			}

			String h = StringUtils.join(header, "\t");
			out.write(h);
			out.newLine();

			// write body
			for (Gene gene: orderedGeneList) {
				String geneName = gene.getName();
                List<String> line = new ArrayList<>(cellBarcodes.size()+1);
                line.add(geneName);
                double [] expressionByGene = dge.getExpression(geneName);
                for (double exp: expressionByGene)
					if (exp==0) // 0 be zero.
                		line.add(Integer.toString(0));
                	else
									line.add(EXPRESSION_FORMATTER.format(exp));

                String b = StringUtils.join(line, "\t");
                out.write(b); out.newLine();
            }
			out.close();
		} catch (IOException e) {
			throw new RuntimeException("Trouble writing " + output.getAbsolutePath(), e);
		}
	}

	/**
	 * Writes the gene expression in a format consumed by tensorQTL.
	 * NOTE: normalize_expression.py reorders the genes by TSS. This is not done here.
	 */
	private void writeExpressionBed(
			final DGEMatrix dge,
			final List<Gene> orderedGeneList,
			final File output
	) {
		try {
			IOUtil.assertFileIsWritable(output);
			final List<String> cellBarcodes = dge.getCellBarcodes();

			final BufferedWriter out =
					new BufferedWriter(new OutputStreamWriter(IOUtil.openFileForWriting(output)));

			// write header
			final List<String> header = new ArrayList<>(cellBarcodes.size() + 6);
			header.add("#chr");
			header.add("start");
			header.add("end");
			header.add("pid");
			for (final String c : cellBarcodes) {
				// get rid of the batch number :
				final String cc = c.split(":")[1];
				header.add(cc);
			}

			final String h = StringUtils.join(header, "\t");
			out.write(h);
			out.newLine();

			// write body
			for (final Gene gene : orderedGeneList) {
				final String geneName = gene.getName();
				final String contig = gene.getContig();
				final int tss = gene.isNegativeStrand() ? gene.getEnd() : gene.getStart();
				final List<String> line = new ArrayList<>(cellBarcodes.size() + 1);
				line.add(contig);
				// BED is 0-based, half-open. So start is TSS - 1.
				line.add(String.valueOf(tss - 1));
				// tensorQTL expects the end to be the same as the start.
				line.add(String.valueOf(tss));
				line.add(geneName);
				final double[] expressionByGene = dge.getExpression(geneName);
				for (final double exp : expressionByGene) {
					if (exp == 0) // 0 be zero.
					{
						line.add(Integer.toString(0));
					} else {
						line.add(EXPRESSION_FORMATTER.format(exp));
					}
				}

				final String b = StringUtils.join(line, "\t");
				out.write(b);
				out.newLine();
			}
			out.close();
		} catch (final IOException e) {
			throw new RuntimeException("Trouble writing " + output.getAbsolutePath(), e);
		}
	}

	private DGEMatrix readAndFilterMatrix (final File f, final Collection <String> rejectedDonors, final String prefix) {
		DGEMatrix m = DGEMatrix.parseDenseFile(f, prefix);
		if (rejectedDonors!=null && rejectedDonors.size()>0) m.toDenseMatrix();
		if (rejectedDonors!=null & !rejectedDonors.isEmpty()) {
			// annoying bit of a hack since we have a prefix.
			Collection <String> rejectedDonorsThis = addPrefixToDonoNames(rejectedDonors, prefix);
			m.removeCellBarcodes(rejectedDonorsThis);
		}
		return m;
	}

	private Collection <String> addPrefixToDonoNames (final Collection<String> donors, final String prefix) {
		List<String> result = new ArrayList<>(donors.size());
		for (String r: donors) {
			String n = prefix+r;
			result.add(n);
		}
		return (result);
	}
	private void validateMetaCellFiles (final List<File> files) {
		for (File f: files)
			IOUtil.assertFileIsReadable(f);
	}

	/**
	 * Sorts genes by the sequence index of the sequence dictionary (instead of the contig string name), then start/end.
	 * @author nemesh
	 */
	private class GeneComparator implements Comparator<Gene> {

		private final SAMSequenceDictionary dict;

		public GeneComparator (final SAMSequenceDictionary dict) {
			this.dict=dict;
		}

        @Override
        public int compare(final Gene o1, final Gene o2) {
        	int ret = IntervalTagComparator.compare(o1, o2, dict);
            if (ret == 0)
                ret = o1.compareTo(o2);
            return ret;
        }
    }

	/** Stock main method. */
	public static void main(final String[] args) {

		System.exit(new PrepareEqtlExpressionData().instanceMain(args));
	}

}
