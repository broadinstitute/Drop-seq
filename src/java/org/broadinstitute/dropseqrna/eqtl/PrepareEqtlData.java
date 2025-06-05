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

package org.broadinstitute.dropseqrna.eqtl;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.tools.DGEMatrix;
import org.broadinstitute.dropseqrna.cmdline.CustomCommandLineValidationHelper;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.StringUtil;
import org.broadinstitute.dropseqrna.utils.AssertSequenceDictionaryIntersection;
import org.broadinstitute.dropseqrna.utils.VCFUtils;

import com.google.common.collect.Sets;
import com.google.common.collect.Sets.SetView;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.nio.PicardHtsPath;

@CommandLineProgramProperties(
        summary = "Generate data for Matrix EQTL and tensorQTL packages",
        oneLineSummary = "Generate data for Matrix EQTL and tensorQTL packages",
        programGroup = DropSeq.class
)
public class PrepareEqtlData extends CommandLineProgram {

    public static final double DEFAULT_REMOVE_PCT_EXPRESSION = 50d;
    public static final double DEFAULT_FRACTION_SAMPLES_PASSING = 0.9;
    public static final double DEFAULT_HWE_PVALUE = 1e-4;
    public static final double DEFAULT_MAF = 0.05;
    public static final int DEFAULT_GQ_THRESHOLD = 30;
    private final Log log = Log.getInstance(PrepareEqtlData.class);

	/**
	 * INPUTS
	 */

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input meta cell expression files.", minElements=1)
	public List<File> META_CELL_FILE;

	@Argument(doc = "The input meta covariate files.", minElements=1)
	public List<File> COVARIATE_FILE;

	@Argument(doc = "Covariate validation stringency.")
	public ValidationStringency COVARIATE_VALIDATION = ValidationStringency.SILENT;

	@Argument(doc = "The input VCF.")
	public PicardHtsPath INPUT_VCF;

	@Argument (doc = "The annotations file that provides gene locations.  Supports GTF and RefFlat format.")
	public File ANNOTATIONS_FILE;

	@Argument (doc="A file with 1 column no header, containing identifiers of donors that should be excluded from the expression output.", optional=true)
	public File REJECTED_DONOR_LIST;

	@Argument(shortName = StandardOptionDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, optional = true,
    doc="Sequences expected to be in ANNOTATIONS_FILE.  If not specified, sequence dictionary is expected to be in same " +
            "directory as ANNOTATIONS_FILE.")
    public File SEQUENCE_DICTIONARY;

	@Argument(
			doc = "If an interval file is supplied, emit SNPs that are included in the interval.  This file is in " +
					"Interval format - tab seperated with fields: chr start end strand name",
			optional = true
	)
	public File INTERVAL_FILE;

	/**
	 * OUTPUTS
	 */
	@Argument(
			doc = "The output expression matrix.  First entry is the gene symbol, columns following that "
					+ "are each donor's expression.  One gene per row.",
			optional = true
	)
	public File EXPRESSION_FILE;

	@Argument(doc = "The output locations of the genes in the expression matrix", optional = true)
	public File GENE_LOCATION_FILE;

	@Argument (doc="The output covariate matrix for the experiment")
	public File COVARIATE_MATRIX;

	@Argument(
			doc="The output genotype matrix. "
					+ "See Matrix eQTL for more details: "
					+ "https://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html",
			optional = true
	)
	public File GENOTYPE_MATRIX;

	@Argument(
			doc="The output snp locations.  "
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
					+ "See Plink 1.9 for more details: "
					+ "https://www.cog-genomics.org/plink/1.9/formats#tped",
			optional = true
	)
	public File PLINK_TPED;

	@Argument(
			doc="The output phenotype expression BED. "
					+ "See tensorQTL for more details: "
					+ "https://github.com/broadinstitute/tensorqtl/tree/master#input-formats",
			optional = true
	)
	public File EXPRESSION_BED_FILE;

	@Argument(doc = "A list of donors which may contain repeats.  " +
			"All donors from the META_CELL_FILES, with the rejected donors from the REJECTED_DONOR_LIST removed.",
			optional = true
	)
	public File OUT_DONOR_LIST;

	/**
	 * OPTIONS
	 */
	@Argument (doc="Should a covariate be added that encodes each batch?  This adds the same covariate parameter to each covariate file, with a string value specific to each input file, "
			+ "representing the experimental batches.")
	public boolean ADD_BATCH_COVARIATE=false;

	@Argument (doc="Each meta cell will have this many total transcripts distributed across all genes.  Set this value to -1 to not normalize donor expression to be"
			+ "new expression = expression gene/sum(all genes) *TRANSCRIPTS_PER_CELL ")
	public Integer TRANSCRIPTS_PER_CELL=100000;

	@Argument (doc="Remove the bottom PCT of expressed genes from the data set. If unset no data is removed.  This calculates the median expression across for each gene "
			+ "across meta cells, then orders genes by expression and calculates percentiles.  The bottom <REMOVE_PCT_EXPRESSION> expressed genes are filtered.", optional=true)
	public Double REMOVE_PCT_EXPRESSION = DEFAULT_REMOVE_PCT_EXPRESSION;

	@Argument(doc= "The minimum genotype quality for a variant.")
	public Integer GQ_THRESHOLD= DEFAULT_GQ_THRESHOLD;

	@Argument (doc="At least <FRACTION_SAMPLES_PASSING> samples must have genotype scores >= GQ_THRESHOLD for the variant in the VCF to be included in the analysis.")
	public double FRACTION_SAMPLES_PASSING = DEFAULT_FRACTION_SAMPLES_PASSING;

	@Argument(doc="HWE pvalue filter.  Calculates HWE based on all donors in VCF")
	public double HWE_PVALUE= DEFAULT_HWE_PVALUE;

	@Argument(doc="Minor allele freqeuncy filter.  Calculates minor allele frequency based on the input set of donors.")
	public double MAF= DEFAULT_MAF;

	@Argument(doc="Filter both very low MAFs and very high MAFs.  Filters variants with MAF < threshold and MAF > (1-threshold)")
	public boolean SYMMETRIC_AF=true;

	@Argument (doc="A list of chromosomes to omit from the analysis.  The default is to omit the Y chromosome and MT.")
	public List<String> IGNORED_CHROMOSOMES= new ArrayList<>(Arrays.asList("Y", "MT"));

	// Arguments below are passed through to PrepareEqtlSnpGeneMap
	@Argument(doc="If set, output file of {SNP, gene} pairs to which eQTL discovery will be restricted. " +
			"If set, at least one of CIS_DIST, NS_INTERVAL_LIST or GS_INTERVAL_LIST must be set.",
			optional = true)
	public File SNP_GENE_MAP;

	@Argument(doc="Consider a SNP this close to a gene.  If not set, do not consider SNPs by distance only.",
			optional = true)
	public Integer CIS_DIST;

	@Argument(doc="List of non-gene-specific intervals to find SNPs to be tested with all genes within NS_CIS_DIST of the SNP.",
			optional = true)
	public File NS_INTERVAL_LIST;

	@Argument(doc="Test genes this close to SNPs in intervals in NS_INTERVAL_LIST.")
	public int NS_CIS_DIST = Integer.MAX_VALUE;

	@Argument(doc="List of gene-specific intervals to find SNPs to be tested with the gene within GS_CIS_DIST of the SNP. " +
			"The name of the interval is the gene for which SNPs in the interval are to be tested.",
			optional = true)
	public File GS_INTERVAL_LIST;

	@Argument(doc="Test genes this close to SNPs in intervals in GS_INTERVAL_LIST.")
	public int GS_CIS_DIST = Integer.MAX_VALUE;


	@Override
	protected String[] customCommandLineValidation() {
		String[] validationErrors = null;
		if (SNP_GENE_MAP != null) {
			validationErrors = makePrepareEqtlSnpGeneMap().customCommandLineValidation();
		}
		return CustomCommandLineValidationHelper.makeValue(super.customCommandLineValidation(), validationErrors);
	}

	private PrepareEqtlSnpGeneMap makePrepareEqtlSnpGeneMap() {
		return new PrepareEqtlSnpGeneMap(
				SNP_GENE_MAP, SNP_LOCATIONS, GENE_LOCATION_FILE, CIS_DIST, NS_INTERVAL_LIST, NS_CIS_DIST,
				GS_INTERVAL_LIST, GS_CIS_DIST, getCommandLine());
	}

	@Override
	public int doWork() {
		// validate inputs/outputs
		validateInputFiles(this.META_CELL_FILE);
		validateInputFiles(this.COVARIATE_FILE);
		IOUtil.assertFileIsReadable(INPUT_VCF.toPath());
		IOUtil.assertFileIsReadable(ANNOTATIONS_FILE);
		if (REJECTED_DONOR_LIST!=null) IOUtil.assertFileIsReadable(REJECTED_DONOR_LIST);

		if (this.EXPRESSION_FILE != null) {
			IOUtil.assertFileIsWritable(this.EXPRESSION_FILE);
			IOUtil.assertFileIsWritable(this.GENE_LOCATION_FILE);
		}
		if (this.EXPRESSION_BED_FILE != null) {
			IOUtil.assertFileIsWritable(this.EXPRESSION_BED_FILE);
		}
		IOUtil.assertFileIsWritable(this.COVARIATE_MATRIX);
		PrepareEqtlGenotypeData.assertGenotypesAreWritable(
				this.GENOTYPE_MATRIX,
				this.SNP_LOCATIONS,
				this.GENOTYPE_BED,
				this.PLINK_TPED
		);

		// validate more inputs as possible.
		if (this.META_CELL_FILE.size()!=this.COVARIATE_FILE.size()) {
			log.error ("Must have the same number of meta cell files and covariate files.");
			System.exit(1);
		}

		// validate that there's an index for the VCF file.
		// otherwise exit.
		if (!VCFUtils.hasIndex(this.INPUT_VCF.toPath())) return 1;

		this.SEQUENCE_DICTIONARY=getSequenceDictionaryFile(this.SEQUENCE_DICTIONARY, this.ANNOTATIONS_FILE);

        AssertSequenceDictionaryIntersection.assertIntersectionVcfBam(
				INPUT_VCF.toPath(), SEQUENCE_DICTIONARY.toPath(), log
		);

        boolean passValidation = validateMetaCellCovariateDonorList(this.META_CELL_FILE, this.COVARIATE_FILE);
        if (!passValidation) {
            log.error("Covariate and meta cell files did not pass validation, exiting!");
            return (1);
        }

        // steps:
		// 1 - prepare expression data
        // 2 - prepare covariate data
		// 3- prepare genotype data
        final File sampleFile = prepareExpression();
		prepareCovariateData(sampleFile);
		prepareGenotypeData(sampleFile);
		prepareSnpMap();
		return 0;
	}


	public static File getSequenceDictionaryFile (final File sequenceDictionaryFile, final File annotationsFile) {
		File result;
		if (sequenceDictionaryFile != null) {
			result = sequenceDictionaryFile;
		    if (!sequenceDictionaryFile.exists()) {
		        throw new RuntimeException(String.format("SEQUENCE_DICTIONARY %s does not exist",
		        		sequenceDictionaryFile.getAbsolutePath()));
            }
        } else {
            String basename = annotationsFile.getName();
		    basename = StringUtil.maybeRemoveSuffix(basename, ".gz");
		    if (basename.endsWith(".gtf")) {
		        basename = StringUtil.maybeRemoveSuffix(basename, ".gtf");
            } else if (basename.endsWith(".refFlat")) {
                basename = StringUtil.maybeRemoveSuffix(basename, ".refFlat");
            }
		    result = new File(annotationsFile.getParentFile(), basename + ".dict");
		    if (!result.exists()) {
		        throw new RuntimeException(String.format("Could not find sequence dictionary for ANNOTATIONS_FILE %s",
		        		annotationsFile.getAbsolutePath()));
            }
        }
		return (result);
	}


	private boolean validateMetaCellCovariateDonorList (final List<File> metaCellFiles, final List<File> covariateFiles) {
		// don't test if there are no covariates.
		if (covariateFiles==null || covariateFiles.isEmpty()) return (true);
		log.info("Validating covariate files and metaCell files have the same donor lists.");

		// don't filter donors.
		boolean pass=true;
		Set<String> donors = Collections.emptySet();
		for (int i=0; i<metaCellFiles.size(); i++) {
			EqtlCovariate cov = EqtlCovariate.parseFile(covariateFiles.get(i), donors, COVARIATE_VALIDATION);
			DGEMatrix m = DGEMatrix.parseDenseFile(metaCellFiles.get(i), null, false);
			Set<String> covDonors = new HashSet<>(cov.donorNames());
			Set<String> metaDonors = new HashSet<>(m.getCellBarcodes());
			if (!covDonors.containsAll(metaDonors) || !metaDonors.containsAll(covDonors)) {
				pass=false;
				log.error("Covariate donors and meta cell donors don't match up for files " + covariateFiles.get(i).getAbsolutePath() + " " + metaCellFiles.get(i).getAbsolutePath());
				SetView<String> difference = Sets.symmetricDifference(covDonors, metaDonors);
				if (!difference.isEmpty())
					log.error("Donors not in both sets:" + difference);
			}
		}
		if (pass) log.info("All files passed validation!");
		return (pass);
	}

	private void prepareGenotypeData (final File sampleFile) {
		// public int processData (final File inputVCF, final File sampleFile, final File genotypeMatrix, final File snpLocations, final double maf, final double hwePvalue, final int gqThreshold,
		// final double fractionDonorsPassing, final List<String> ignoredChromosomes, final boolean symmetricAlleleFreq) {
		PrepareEqtlGenotypeData prepGenotypeData = new PrepareEqtlGenotypeData();
		prepGenotypeData.VALIDATION_STRINGENCY = VALIDATION_STRINGENCY;
		prepGenotypeData.processData(this.INPUT_VCF, this.INTERVAL_FILE, sampleFile, this.GENOTYPE_MATRIX,
				this.SNP_LOCATIONS, this.GENOTYPE_BED, this.PLINK_TPED, this.MAF, this.HWE_PVALUE, this.GQ_THRESHOLD,
				this.FRACTION_SAMPLES_PASSING, this.IGNORED_CHROMOSOMES, this.SYMMETRIC_AF);
	}

	private void prepareSnpMap() {
		if (SNP_GENE_MAP != null) {
			makePrepareEqtlSnpGeneMap().processData();
		}
	}

	private File prepareExpression() {
		File sampleFile = OUT_DONOR_LIST;
		try {
			if (sampleFile == null) {
				sampleFile = File.createTempFile("tmp_donor_list", ".txt", TMP_DIR.getFirst());
			}
		} catch (IOException e) {
			throw new RuntimeException("Could not create temp file for donor list.", e);
		}
		PrepareEqtlExpressionData prepExpression = new PrepareEqtlExpressionData();
		prepExpression.VALIDATION_STRINGENCY = VALIDATION_STRINGENCY;
		prepExpression.processData(
				this.META_CELL_FILE, this.REJECTED_DONOR_LIST, this.ANNOTATIONS_FILE, SEQUENCE_DICTIONARY,
				this.EXPRESSION_FILE, this.GENE_LOCATION_FILE, this.EXPRESSION_BED_FILE, sampleFile,
				this.TRANSCRIPTS_PER_CELL, this.IGNORED_CHROMOSOMES, this.REMOVE_PCT_EXPRESSION
		);
		return sampleFile;
	}

	private void prepareCovariateData (final File sampleFile) {
		List<String> requestedDonorList =ParseBarcodeFile.readCellBarcodeFile(sampleFile);
		PrepareEqtlCovariates prepCovariates = new PrepareEqtlCovariates();
		prepCovariates.VALIDATION_STRINGENCY = VALIDATION_STRINGENCY;
		EqtlCovariate result = prepCovariates.mergeCovariates(this.COVARIATE_FILE, this.REJECTED_DONOR_LIST, this.ADD_BATCH_COVARIATE);
		result.writeFile(this.COVARIATE_MATRIX, requestedDonorList);
	}

	private void validateInputFiles (final List<File> files) {
		for (File f: files)
			IOUtil.assertFileIsReadable(f);
	}

	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new PrepareEqtlData().instanceMain(args));
	}
}
