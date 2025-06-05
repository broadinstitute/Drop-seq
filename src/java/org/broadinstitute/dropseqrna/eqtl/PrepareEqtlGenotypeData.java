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

import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import htsjdk.tribble.annotation.Strand;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.annotation.GeneAnnotationReader;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.IteratorOfIterators;
import org.broadinstitute.dropseqrna.utils.VariantContextProgressLoggerIterator;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintWriter;
import org.broadinstitute.dropseqrna.vcftools.SampleAssignmentVCFUtils;
import org.broadinstitute.dropseqrna.vcftools.filters.*;
import picard.annotation.Gene;
import picard.cmdline.CommandLineProgram;
import picard.nio.PicardHtsPath;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.PrintWriter;
import java.util.*;

@CommandLineProgramProperties(
        summary = "Generate genotype matrix data for Matrix EQTL and tensorQTL packages",
        oneLineSummary = "Generate genotype matrix data for Matrix EQTL and tensorQTL packages",
        programGroup = DropSeq.class
)
public class PrepareEqtlGenotypeData extends CommandLineProgram {

	private final Log log = Log.getInstance(PrepareEqtlGenotypeData.class);

	@Argument(doc = "The input VCF.")
	public PicardHtsPath INPUT_VCF;

	@Argument(
			doc = "If an interval file is supplied, emit SNPs that are included in the interval.  This file is in " +
					"Interval format - tab seperated with fields: chr start end strand name",
			optional = true
	)
	public File INTERVAL_FILE;

	@Argument(doc = "A file containing a list of samples to filter the VCF by.  Has 1 column, 1 entry per row.  Each entry is a single sample name from the VCF.  Samples can be repeated, yielding repeated columns in the GENOTYPE_MATRIX output")
	public File SAMPLE_FILE;

	@Argument(doc= "The minimum genotype quality for a variant.")
	public Integer GQ_THRESHOLD=30;

	@Argument (doc="At least <FRACTION_SAMPLES_PASSING> samples must have genotype scores >= GQ_THRESHOLD for the variant in the VCF to be included in the analysis.")
	public double FRACTION_SAMPLES_PASSING=0.9;

	@Argument(doc="HWE pvalue filter.  Calculates HWE based on all donors in VCF")
	public double HWE_PVALUE=1e-4;

	@Argument(doc="Minor allele freqeuncy filter.  Calculates minor allele frequency based on the input set of donors.")
	public double MAF=0.05;

	@Argument(doc="Filter both very low MAFs and very high MAFs.  Filters variants with MAF < threshold and MAF > (1-threshold)")
	public boolean SYMMETRIC_AF=true;

	@Argument(
			doc="The output genotype matrix. "
					+ "At least one of GENOTYPE_MATRIX, GENOTYPE_BED, and PLINK_TPED must be set. "
					+ "When GENOTYPE_MATRIX is set SNP_LOCATIONS must be set. "
					+ "See Matrix eQTL for more details: "
					+ "https://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html",
			optional = true
	)
	public File GENOTYPE_MATRIX;

	@Argument(
			doc="The output snp locations. "
					+ "Must be set when GENOTYPE_MATRIX is set. "
					+ "See Matrix eQTL for more details: "
					+ "https://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html",
			optional = true
	)
	public File SNP_LOCATIONS;

	@Argument(
			doc = "The output tensorQTL genotype bed. "
					+ "At least one of GENOTYPE_MATRIX, GENOTYPE_BED, and PLINK_TPED must be set. "
					+ "See this commit for more details: "
					+ "https://github.com/broadinstitute/tensorqtl/commit/94a7553",
			optional = true
	)
	public File GENOTYPE_BED;

	@Argument(
			doc="The output plink tped. "
					+ "At least one of GENOTYPE_MATRIX, GENOTYPE_BED, and PLINK_TPED must be set. "
					+ "See Plink 1.9 for more details: "
					+ "https://www.cog-genomics.org/plink/1.9/formats#tped",
			optional = true
	)
	public File PLINK_TPED;

	@Argument (doc="A list of chromosomes to omit from the analysis.  The default is to omit the Y chromosome and MT.")
	public List<String> IGNORED_CHROMOSOMES= new ArrayList<>(Arrays.asList("Y", "MT"));

	@Argument (doc="Log every rejected variant to the console with the reason it was rejected.", optional=true)
	public boolean VERBOSE=false;

	/**
	 * File extension for .tped files.
	 *
	 * See PLINK 1.9 for more details:
	 * <a href="https://www.cog-genomics.org/plink/1.9/formats#tped">https://www.cog-genomics.org/plink/1.9/formats#tped</a>
	 */
	public static final String TPED_FILE_EXTENSION = ".tped";

	/**
	 * File extension for .tfam files.
	 *
	 * See PLINK 1.9 for more details:
	 * <a href="https://www.cog-genomics.org/plink/1.9/formats#tfam">https://www.cog-genomics.org/plink/1.9/formats#tfam</a>
	 */
	public static final String TFAM_FILE_EXTENSION = ".tfam";

	/**
	 * File extension for reference a2-allele files.
	 *
	 * See PLINK 1.9 for more details:
	 * <a href="https://www.cog-genomics.org/plink/1.9/data#ax_allele">https://www.cog-genomics.org/plink/1.9/data#ax_allele</a>
	 */
	public static final String REF_ALLELE_FILE_EXTENSION = ".ref_allele";

	/**
	 * A "no call" entry for a PED allele.
	 *
	 * Via: "'0' = no call" in
	 * <a href="https://www.cog-genomics.org/plink/1.9/formats#ped">https://www.cog-genomics.org/plink/1.9/formats#ped</a>.
	 */
	private static final String PLINK_NO_CALL = "0";

	private static final int IO_BUFFER_SIZE = 10 * 1024 * 1024;

	private static final String MATRIX_EQTL_FIELD_SEPARATOR = "\t";

	private static final String TENSORQTL_FIELD_SEPARATOR = "\t";

	private static final String PLINK_FIELD_SEPARATOR = " ";

	/**
	 *
	 * The end to end processing of this class encapsulated in a single method to allow pipelining.
	 * See program options for details.
	 */
	@Override
	public int doWork() {
		return processData(
				this.INPUT_VCF, this.INTERVAL_FILE, this.SAMPLE_FILE, this.GENOTYPE_MATRIX, this.SNP_LOCATIONS,
				this.GENOTYPE_BED, this.PLINK_TPED, this.MAF, this.HWE_PVALUE, this.GQ_THRESHOLD,
				this.FRACTION_SAMPLES_PASSING, this.IGNORED_CHROMOSOMES, this.SYMMETRIC_AF);
	}

	/**
	 * /**
	 * The end to end processing of this class encapsulated in a single method to allow pipelining.
	 * See program options for details.
	 */
	public int processData(
			final PicardHtsPath inputVCF, final File intervalFile, final File sampleFile, final File genotypeMatrix,
			final File snpLocations, final File genotypeBed, final File plinkTped, final double maf,
			final double hwePvalue, final int gqThreshold, final double fractionDonorsPassing,
			final List<String> ignoredChromosomes, final boolean symmetricAlleleFreq) {

		IOUtil.assertFileIsReadable(inputVCF.toPath());
		IOUtil.assertFileIsReadable(sampleFile);
		PrepareEqtlGenotypeData.assertGenotypesAreWritable(
				genotypeMatrix,
				snpLocations,
				genotypeBed,
				plinkTped
		);

		// due to a weird bug in picard, I think you actually need to have an index...
		final VCFFileReader vcfReader = new VCFFileReader(inputVCF.toPath(), true);

		// you need all the requested donors present in the VCF file, or something has gone wrong.
		List<String> requestedDonorList =ParseBarcodeFile.readCellBarcodeFile(sampleFile);
		if (!validateAllDonorsInVCF(vcfReader, requestedDonorList))
			return 1;

		PeekableIterator<VariantContext> iter = getVCFIterator (vcfReader, intervalFile, requestedDonorList, maf, hwePvalue, gqThreshold,
				fractionDonorsPassing, ignoredChromosomes, symmetricAlleleFreq);

		PrintWriter snpLocationWriter = null;
		PrintWriter genotypeMatrixWriter = null;
		if (genotypeMatrix != null) {
			snpLocationWriter = bufferedWriter(snpLocations);
			genotypeMatrixWriter = bufferedWriter(genotypeMatrix);

			// prep headers
			writeSNPLocationHeader(snpLocationWriter);
			writeGenotypeMatrixHeader(genotypeMatrixWriter, requestedDonorList);
		}

		PrintWriter genotypeBedWriter = null;
		if (genotypeBed != null) {
			genotypeBedWriter = bufferedWriter(genotypeBed);

			writeGenotypeBedHeader(genotypeBedWriter, requestedDonorList);
		}

		PrintWriter plinkTpedWriter = null;
		PrintWriter plinkRefAlleleWriter = null;
		PrintWriter plinkTfamWriter = null;
		if (plinkTped != null) {
			plinkTpedWriter = bufferedWriter(plinkTped);
			plinkRefAlleleWriter =
					bufferedWriter(resolvePlinkSibling(plinkTped, REF_ALLELE_FILE_EXTENSION));
			plinkTfamWriter = bufferedWriter(resolvePlinkSibling(plinkTped, TFAM_FILE_EXTENSION));

			writePlinkTfam(plinkTfamWriter, requestedDonorList);
		}

			while (iter.hasNext()) {
				VariantContext vc = iter.next();
				writeVariant(
						vc, requestedDonorList, gqThreshold, snpLocationWriter, genotypeMatrixWriter,
						genotypeBedWriter, plinkTpedWriter, plinkRefAlleleWriter
				);
			}

		CloserUtil.close(snpLocationWriter);
		CloserUtil.close(genotypeMatrixWriter);
		CloserUtil.close(genotypeBedWriter);
		CloserUtil.close(plinkTpedWriter);
		CloserUtil.close(plinkRefAlleleWriter);
		CloserUtil.close(plinkTfamWriter);
		CloserUtil.close(vcfReader);
		return 0;
	}

	private static PrintWriter bufferedWriter(final File output) {
		return new ErrorCheckingPrintWriter(
				new BufferedOutputStream(
						IOUtil.openFileForWriting(output),
						IO_BUFFER_SIZE
				)
		);
	}

	/**
	 * Convert data to a matrix of dosages of the alternate allele
	 */
	private void writeVariant(final VariantContext vc, final List<String> requestedDonorList,
			final int gqThreshold, final PrintWriter snpLocationWriter,
							  final PrintWriter genotypeMatrixWriter, final PrintWriter genotypeBedWriter,
							  final PrintWriter plinkTpedWriter, final PrintWriter plinkRefAlleleWriter) {
		Allele refAllele = vc.getReference();
		Allele altAllele = vc.getAltAlleleWithHighestAlleleCount();
		final String variantContig = vc.getContig();
		final String variantStart = String.valueOf(vc.getStart());
		// need to be able to encode the variant name if the alt allele is null.
		String variantName = variantContig + ":" + variantStart + ":" + refAllele.getBaseString() + ":";
		if (altAllele!=null)
			variantName+=altAllele.getBaseString();
		else
			variantName+="-";

		final ArrayList<String> dosAlts = new ArrayList<>();
		for (String donor: requestedDonorList) {
			Genotype g = vc.getGenotype(donor);
			// low quality / not called genotype.
			if (g.getGQ() < gqThreshold | !g.isCalled()) {
				dosAlts.add("NA");
			} else {
				int dosAlt = g.countAllele(altAllele);
				dosAlts.add(Integer.toString(dosAlt));
			}
		}

		if (genotypeMatrixWriter != null) {
			// write out the line to the genotype matrix.
			// write out the line to the locations
			writeGenotypeMatrixLine(genotypeMatrixWriter, snpLocationWriter, variantName, vc, dosAlts);
		}

		if (genotypeBedWriter != null) {
			// write out the line to the genotype bed.
			writeGenotypeBedLine(genotypeBedWriter, variantName, vc, dosAlts);
		}

		if (plinkTpedWriter != null) {
			// Write the reference alleles
			final String refBase = refAllele.isNoCall() ? PLINK_NO_CALL : refAllele.getBaseString();
			final String altBase =
					altAllele == null || altAllele.isNoCall() ? PLINK_NO_CALL : altAllele.getBaseString();
			plinkRefAlleleWriter.printf("%s\t%s%n", variantName, refBase);
			writePlinkTpedLine(plinkTpedWriter, variantName, vc, refBase, altBase, dosAlts);
		}
	}

	private void writeGenotypeMatrixHeader(
			final PrintWriter genotypeMatrixWriter,
			final List<String> donors) {
		List<String> header = new ArrayList<>();
		header.add("id");
		header.addAll(donors);
		genotypeMatrixWriter.println(StringUtils.join(header, MATRIX_EQTL_FIELD_SEPARATOR));
	}

	private void writeSNPLocationHeader(final PrintWriter snpLocationWriter) {
		String [] header ={"snp", "chr", "pos", "id", "end"};
		snpLocationWriter.println(StringUtils.join(header, MATRIX_EQTL_FIELD_SEPARATOR));
	}

	private void writeGenotypeBedHeader(
			final PrintWriter genotypeMatrixWriter,
			final List<String> donors) {
		List<String> header = new ArrayList<>();
		header.add("#chr");
		header.add("start");
		header.add("end");
		header.add("pid");
		header.addAll(donors);
		genotypeMatrixWriter.println(StringUtils.join(header, TENSORQTL_FIELD_SEPARATOR));
	}

	private static void writePlinkTfam(final PrintWriter plinkTfamWriter, final List<String> donors) {
		for (final String donor : donors) {
      /*
      .tfam is the same as .fam, but it accompanies .tped instead of .ped.

      Columns via:
       - https://www.cog-genomics.org/plink/1.9/formats#fam

       1. Family ID: 'FID'
       2. Within-family ID: 'IID'
       3. Within-family ID of father: '0' if father isn't in dataset
       4. Within-family ID of mother: '0' if mother isn't in dataset
       5. Sex code: '0' = unknown
       6. Phenotype value: '-9' = missing data if case/control
       */
			final String[] line = {"0", donor, "0", "0", "0", "-9"};
			plinkTfamWriter.println(StringUtils.join(line, PLINK_FIELD_SEPARATOR));
		}
	}

	private void writeGenotypeMatrixLine(
			final PrintWriter genotypeMatrixWriter,
			final PrintWriter snpLocationWriter,
			final String variantName,
			final VariantContext vc,
			final List<String> dosAlts) {
		final List<String> genotypeMatrixLine = new ArrayList<>();
		genotypeMatrixLine.add(variantName);
		genotypeMatrixLine.addAll(dosAlts);
		genotypeMatrixWriter.println(StringUtils.join(genotypeMatrixLine, MATRIX_EQTL_FIELD_SEPARATOR));

		// write out the line to the locations
		final String[] snpLoc = {
				variantName,
				vc.getContig(),
				Integer.toString(vc.getStart()), vc.getID(),
				Integer.toString(vc.getEnd())
		};
		snpLocationWriter.println(StringUtils.join(snpLoc, MATRIX_EQTL_FIELD_SEPARATOR));
	}

	private void writeGenotypeBedLine(
			final PrintWriter genotypeBedWriter,
			final String variantName,
			final VariantContext vc,
			final List<String> dosAlts) {
		final int start = vc.getStart();
		final List<String> genotypeBedLine = new ArrayList<>();
		genotypeBedLine.add(vc.getContig());
		// BED is 0-based, half-open. So start - 1.
		genotypeBedLine.add(Integer.toString(start - 1));
		// tensorQTL expects the end to be the same as the start.
		genotypeBedLine.add(Integer.toString(start));
		genotypeBedLine.add(variantName);
		genotypeBedLine.addAll(dosAlts);
		genotypeBedWriter.println(StringUtils.join(genotypeBedLine, TENSORQTL_FIELD_SEPARATOR));
	}

	private void writePlinkTpedLine(
			final PrintWriter plinkTpedWriter,
			final String variantName,
			final VariantContext vc,
			final String refBase,
			final String altBase,
			final List<String> dosAlts) {
		/*
		First four columns of .tped are the same as .map.

		Columns via:
		- https://www.cog-genomics.org/plink/1.9/formats#tped
		- https://www.cog-genomics.org/plink/1.9/formats#map

		1. Chromosome code
		2. Variant identifier
		3. Position in morgans or centimorgans (optional; also safe to use dummy value of '0')
		4. Base-pair coordinate
		*/
		final List<String> plinkTpedLine = new ArrayList<>();
		plinkTpedLine.add(vc.getContig());
		plinkTpedLine.add(variantName);
		plinkTpedLine.add("0");
		plinkTpedLine.add(Integer.toString(vc.getStart()));

		for (final String dosAlt : dosAlts) {
			switch (dosAlt) {
				case "NA":
					plinkTpedLine.add(PLINK_NO_CALL);
					plinkTpedLine.add(PLINK_NO_CALL);
					break;
				case "0":
					plinkTpedLine.add(refBase);
					plinkTpedLine.add(refBase);
					break;
				case "1":
					plinkTpedLine.add(altBase);
					plinkTpedLine.add(refBase);
					break;
				case "2":
					plinkTpedLine.add(altBase);
					plinkTpedLine.add(altBase);
					break;
				default:
					final String message = String.format(
							"Unexpected alt allele count %s: %s %s", dosAlt, variantName, dosAlt);
					switch (VALIDATION_STRINGENCY) {
						case STRICT: throw new SAMException(message);
						case LENIENT: log.warn(message + " -- using two alt alleles"); break;
						case SILENT: break;
					}
					plinkTpedLine.add(altBase);
					plinkTpedLine.add(altBase);
					break;
			}
		}

		plinkTpedWriter.println(StringUtils.join(plinkTpedLine, PLINK_FIELD_SEPARATOR));
	}

	/**
	 * All donors requested must be contained in the VCF.  The donor list may have repeated donors.
	 */
	private boolean validateAllDonorsInVCF (final VCFFileReader vcfReader, final List<String> requestedDonorList) {
		Set<String> samplesInVCF = new HashSet<>(vcfReader.getFileHeader().getSampleNameToOffset().keySet());
		Set<String> requestedDonorSet = new HashSet<>(requestedDonorList);
		// are there requested donors not in the VCF?
		requestedDonorSet.removeAll(samplesInVCF);
		// if there are, log them and error out.
		if (!requestedDonorSet.isEmpty()) {
			log.error("Not all requested donors found in VCF!  Program will exit.");
			for (String s: requestedDonorSet)
				log.error("Missing Donor: " + s);
			return false;
		}
		return true;

	}

	PeekableIterator<VariantContext> getVCFIterator(
			final VCFFileReader vcfReader, final File intervalFile, final List<String> requestedDonorList,
			final double mafThreshold, final double hwePValueThreshold, final int gqThreshold,
			final double fracSamplesPassing, final List<String> ignoredChromosomes, final boolean symmetric
	) {
		log.info("Filtering variants with a MAF of less than [" + mafThreshold+ "]");
		log.info("Filtering variants with HWE pvalue less than [" + hwePValueThreshold +"]" );

		Iterator<VariantContext> iter;

		if (intervalFile == null) {
			iter = vcfReader.iterator();
		} else {
			iter = IteratorOfIterators.fromIterator(
					IntervalList.fromFile(intervalFile).iterator(),
					interval -> vcfReader.query(interval.getContig(), interval.getStart(), interval.getEnd())
			);
		}

		iter = new VariantContextProgressLoggerIterator(iter, new ProgressLogger(log));

		iter = new SampleAssignmentVCFUtils.VariantContextSampleTransformer (iter, new HashSet<>(requestedDonorList));

		// filter on the variant to remove flagged variants (not PASS) with 2 alleles maximum, with unlimited length, drop monomorphic SNPs,
		iter = new MonomorphicVariantContextFilter(iter, requestedDonorList);

		iter = new SimpleDiploidVariantContextFilter(iter, false, true, 2, null, VERBOSE);
		// filter on genotype call rate
		iter = new CallRateVariantContextFilter(iter, gqThreshold, fracSamplesPassing, requestedDonorList, VERBOSE);
		// filter on chromosomes.
		iter = new PeekableIterator<>(new ChromosomeVariantFilter(iter, ignoredChromosomes));

		// need filters on MAF and HWE
		iter = new HardyWeinbergVariantContextFilter(iter, hwePValueThreshold, VERBOSE);
		iter = new MinorAlleleFreqVariantContextFilter(iter, mafThreshold, symmetric, requestedDonorList, gqThreshold, VERBOSE);
		return new PeekableIterator<>(iter);
	}

   /**
    * Get the intervals of genes from an annotation File.
    * If the annotations file is null, return null
    * @param annotationFile The annotation file to parse.
    * @return An IntervalList, one element per gene defining the bounds of the genes in the annotation File.
    */
	public static IntervalList parseGTFToIntervalList(
			File annotationFile, File sequenceDictionaryFile, final ValidationStringency validationStringency
	) {
		if (annotationFile == null) return null;
		SamReader reader = SamReaderFactory.makeDefault().open(sequenceDictionaryFile);
		SAMSequenceDictionary dict = reader.getFileHeader().getSequenceDictionary();
		CloserUtil.close(reader);

		OverlapDetector<Gene> od = GeneAnnotationReader.loadAnnotationsFile(annotationFile, dict, validationStringency);
		Set<Gene> genes = od.getAll();
		IntervalList intervalList = new IntervalList(dict);
		genes.forEach(x -> intervalList.add(
				new Interval(x.getContig(), x.getStart(), x.getEnd(), x.getStrand() == Strand.REVERSE, x.getName())));
		return (intervalList);

	}

	public static void assertGenotypesAreWritable(
			final File genotypeMatrix,
			final File snpLocations,
			final File genotypeBed,
			final File plinkTped) {
		if (genotypeMatrix == null && genotypeBed == null && plinkTped == null) {
			throw new SAMException("At least one of genotypeMatrix, genotypeBed, and plinkTped should be specified");
		}
		assertGenotypeMatrixIsWritable(genotypeMatrix, snpLocations);
		assertGenotypeBedIsWritable(genotypeBed);
		assertTpedIsWritable(plinkTped);
	}

	private static void assertGenotypeMatrixIsWritable(
			final File genotypeMatrix,
			final File snpLocations) {
		if (genotypeMatrix == null && snpLocations != null) {
			throw new SAMException("If snpLocations is supplied, must also supply genotypeMatrix");
		}

		if (genotypeMatrix != null) {
			if (snpLocations == null) {
				throw new SAMException("If genotypeMatrix is supplied, must also supply snpLocations");
			}
			IOUtil.assertFileIsWritable(genotypeMatrix);
			IOUtil.assertFileIsWritable(snpLocations);
		}
	}

	private static void assertGenotypeBedIsWritable(
			final File genotypeBed) {
		if (genotypeBed != null) {
			IOUtil.assertFileIsWritable(genotypeBed);
		}
	}

	private static void assertTpedIsWritable(final File plinkTped) {
		if (plinkTped != null) {
			if (!plinkTped.getName().endsWith(TPED_FILE_EXTENSION)) {
				throw new IllegalArgumentException(
						String.format(
								"Output tped file %s does not have the extension %s",
								plinkTped.getAbsolutePath(),
								TPED_FILE_EXTENSION
						)
				);
			}
			final File plinkTfam = resolvePlinkSibling(plinkTped, TFAM_FILE_EXTENSION);
			final File plinkRefAllele = resolvePlinkSibling(plinkTped, REF_ALLELE_FILE_EXTENSION);
			IOUtil.assertFileIsWritable(plinkTfam);
			IOUtil.assertFileIsWritable(plinkRefAllele);
		}
	}

	public static File resolvePlinkSibling(final File plinkFile, final String extension) {
		final String basename = IOUtil.basename(plinkFile);
		return new File(plinkFile.getParentFile(), basename + extension);
	}

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new PrepareEqtlGenotypeData().instanceMain(args));
	}
}
