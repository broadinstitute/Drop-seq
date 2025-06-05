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

import com.google.common.collect.Sets.SetView;

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.tools.DGEMatrix;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;

import org.broadinstitute.dropseqrna.eqtl.EqtlCovariate;
import org.broadinstitute.dropseqrna.eqtl.PrepareEqtlData;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.TabbedInputParser;

import java.io.File;
import java.text.DecimalFormat;
import java.util.*;

@CommandLineProgramProperties(
        summary = "Generate covariate matrix data for Matrix EQTL and tensorQTL packages",
        oneLineSummary = "Generate covariate matrix data for Matrix EQTL and tensorQTL packages",
        programGroup = DropSeq.class
)
public class PrepareEqtlCovariates extends CommandLineProgram {

	private final Log log = Log.getInstance(PrepareEqtlCovariates.class);

	@Argument (doc="An input meta cell expression file.  This option will generate an empty covariate file with appropriate headers, emitted to the MERGED_OUTPUT parameter.", mutex={"COVARIATE_FILES"})
	public File META_CELL_FILE;

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input meta covariate files.", mutex={"META_CELL_FILE"})
	public List<File> COVARIATE_FILES;

	@Argument (doc="A file containing columns of covariates and rows of donors.  The first column contains the donor name, and subsequent columns the covariate values.  "
			+ "Each row contains the info for a single donor.  Has a header row containing the donor ID coloumn + covariate names.  "
			+ "If supplied, adds covariates to the output covariate matrix.  This can be specified multiple times.", optional=true)
	public List<File> COVARIATE_REFERENCE_FILE;

	@Argument(doc = "Covariate validation stringency.")
	public ValidationStringency COVARIATE_VALIDATION = ValidationStringency.SILENT;

	@Argument (doc="Add the covariate(s) from the COVARIATE_REFERENCE_FILE to the output covariate matrix.  "
			+ "Can be set multiple times to add multiple columns.", optional=true, mutex=("REFERENCE_COVARIATE_LIST_FILE"))
	public List<String> REFERENCE_COVARIATE;

	@Argument(doc="Use the covariate names in this file that intersect the COVARIATE_REFERENCE_FILE name as covariates in the output.  This takes the place of calling the "
			+ "REFERENCE_COVARIATE argument multiple times.  This can also interpret reserved covariate words [FRACTION_X | library_size | batch_effect]." , mutex=("REFERENCE_COVARIATE"))
	public File REFERENCE_COVARIATE_LIST_FILE;
	
	@Argument (doc="The output covariate matrix for the experiment", optional=false)
	public File OUTPUT;

	@Argument (doc="Should a covariate be added that encodes each batch?  Only used in combination with COVARIATE_FILES parameter.")
	public boolean ADD_BATCH_COVARIATE=false;

	@Argument (doc="Adds a library size covariate that encodes the number of total UMIs for the donor.  Only used when a generating a new covariates file where the META_CELL_FILE is the input.")
	public boolean ADD_LIBRARY_SIZE_COVARIATE=false;
	
	@Argument (doc="Add a covariate for X reactivation.  If this is true, then the following parameters must be set: ANNOTATIONS_FILE, CONTIG_GROUP_FILE")
	public boolean ADD_X_REACTIVATION_COVARIATE=false;
	
	@Argument (doc = "The annotations file that provides gene locations.  Supports GTF and RefFlat format.", optional=true)
	public File ANNOTATIONS_FILE;
	
	@Argument (doc="A YAML file containing the annotation groups each contig belongs to.  The file has a list of contig names, each of which has a list of annotation groups that contig belongs to.", optional=true)
	public File CONTIG_GROUP_FILE;
	
	@Argument(shortName = StandardOptionDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME, optional = true,
    doc="Sequences expected to be in ANNOTATIONS_FILE.  If not specified, sequence dictionary is expected to be in same " +
		            "directory as ANNOTATIONS_FILE.")
	public File SEQUENCE_DICTIONARY;
	
	@Argument (doc="A list of genes that escape X inactivation.  Optional to provide when calculating X reactivation.", optional=true)
	public File ESCAPE_GENES_FILE=null;

	@Argument (doc="A file with 1 column no header, containing identifiers of donors that should be excluded from the expression output. Only used in combination with COVARIATE_FILES parameter", optional=true)
	public File REJECTED_DONOR_LIST;

	private String BATCH_EFFECT_ATTRIBUTE="batch_effect";
	private String LIBRARY_SIZE_ATTRIBUTE="library_size";

	@Override
	public int doWork() {
		// if there's a file containing one or more covariates to include, then read this in and set the variable.
		processCovarList();
		
		if (ADD_X_REACTIVATION_COVARIATE) {
			if (ANNOTATIONS_FILE==null) log.error("If ADD_X_REACTIVATION_COVARIATE is true, must supply ANNOTATIONS_FILE");
			if (CONTIG_GROUP_FILE==null) log.error("If CONTIG_GROUP_FILE is true, must supply ANNOTATIONS_FILE");
			if (ANNOTATIONS_FILE==null || CONTIG_GROUP_FILE==null) 
				return(1);
			IOUtil.assertFileIsReadable(this.ANNOTATIONS_FILE);
			IOUtil.assertFileIsReadable(this.CONTIG_GROUP_FILE);
			if (ESCAPE_GENES_FILE!=null) IOUtil.assertFileIsReadable(ESCAPE_GENES_FILE);
			this.SEQUENCE_DICTIONARY= PrepareEqtlData.getSequenceDictionaryFile(this.SEQUENCE_DICTIONARY, this.ANNOTATIONS_FILE);
		}
		
		if (REFERENCE_COVARIATE.size()>0 && this.COVARIATE_REFERENCE_FILE==null) {
			log.error("If REFERENCE_COVARIATE parameter is applied, COVARIATE_REFERENCE_FILE can not be null.");
			return 1;
		}
		
		if (this.COVARIATE_REFERENCE_FILE.size()>0) {
			for (File f: this.COVARIATE_REFERENCE_FILE)
				IOUtil.assertFileIsReadable(f);
		}
		
		IOUtil.assertFileIsWritable(this.OUTPUT);
		if (this.REJECTED_DONOR_LIST!=null) IOUtil.assertFileIsReadable(this.REJECTED_DONOR_LIST);
		if (COVARIATE_FILES.size()==0 & META_CELL_FILE==null) {
			log.info("Must supply at least one covariate file or meta cell file.");
			return 1;
		}
		if (!this.COVARIATE_FILES.isEmpty())
			for (File f: this.COVARIATE_FILES)
				IOUtil.assertFileIsReadable(f);

		if (this.META_CELL_FILE!=null) IOUtil.assertFileIsReadable(this.META_CELL_FILE);
		if (this.META_CELL_FILE==null & this.ADD_LIBRARY_SIZE_COVARIATE) {
			log.error("Can only add library size covariate to a new covariate file.");
			return 1;
		}

		if (this.META_CELL_FILE!=null && this.COVARIATE_FILES.isEmpty()) {
			EqtlCovariate result = generateEmptyCovariatesFile(this.META_CELL_FILE);
			// other work.
			if (this.ADD_LIBRARY_SIZE_COVARIATE)
				result=addLibrarySize(result, this.META_CELL_FILE);
			if (this.COVARIATE_REFERENCE_FILE!=null) {
				result=addReferenceCovariates(this.COVARIATE_REFERENCE_FILE, this.REFERENCE_COVARIATE, result);
				if (result==null) return 1; // catch a problem with adding reference covariates and exit.
			if (this.ADD_X_REACTIVATION_COVARIATE)
				result=addXReactivationCovariate(result, this.META_CELL_FILE, this.ANNOTATIONS_FILE, CONTIG_GROUP_FILE, SEQUENCE_DICTIONARY, ESCAPE_GENES_FILE);
			}
			result.writeFile(this.OUTPUT);
			return 0;
		}

		// otherwise, you're merging multiple covariate files together, removing donors, and possibly adding a batch effect.
		EqtlCovariate result = mergeCovariates(this.COVARIATE_FILES, this.REJECTED_DONOR_LIST, this.ADD_BATCH_COVARIATE);
		if (this.COVARIATE_REFERENCE_FILE!=null) {
			result=addReferenceCovariates(this.COVARIATE_REFERENCE_FILE, this.REFERENCE_COVARIATE, result);
			if (result==null) return 1; // catch a problem with adding reference covariates and exit.
		}
		result.writeFile(this.OUTPUT);
		return 0;
	}
	
	private void processCovarList  () {
		if (REFERENCE_COVARIATE_LIST_FILE==null) return;
		
		IOUtil.assertFileIsReadable(REFERENCE_COVARIATE_LIST_FILE);
		List<String> covarList = ParseBarcodeFile.readCellBarcodeFile(REFERENCE_COVARIATE_LIST_FILE);
		if (covarList.contains(this.BATCH_EFFECT_ATTRIBUTE)) {
			this.ADD_BATCH_COVARIATE=true;
			covarList.remove(this.BATCH_EFFECT_ATTRIBUTE);
		}
		if (covarList.contains(this.LIBRARY_SIZE_ATTRIBUTE)) { 
			this.ADD_LIBRARY_SIZE_COVARIATE=true;
			covarList.remove(this.LIBRARY_SIZE_ATTRIBUTE);
		}
		if (covarList.contains(CalculateXReactivationCovariate.FRACTION_X_COVAR_NAME)) {			
			this.ADD_X_REACTIVATION_COVARIATE=true;
			covarList.remove(CalculateXReactivationCovariate.FRACTION_X_COVAR_NAME);
		}
		this.REFERENCE_COVARIATE=covarList;
		
	}
	
	
	public EqtlCovariate addXReactivationCovariate(EqtlCovariate covars, File metaCellFile, File annotationsFile, File contigGroupFile, File sequenceDictionaryFile, File escapeGenesFile) {
		CalculateXReactivationCovariate calc = new CalculateXReactivationCovariate();
		EqtlCovariate result = calc.getXReactivationCovariate(metaCellFile, annotationsFile, contigGroupFile,
				sequenceDictionaryFile, escapeGenesFile, VALIDATION_STRINGENCY);
		covars.mergeAttributes(result);				
		return covars;
	}
	
	public EqtlCovariate addReferenceCovariates(final List<File> covarFiles, final List<String> covars, final EqtlCovariate result) {
		DonorCovariates dc =parse(covarFiles);
		
		List<String> donorNames = result.donorNames();
		Set<String> refDonorNames = new HashSet<>(dc.getDonorNames());
		SetView<String> difference = com.google.common.collect.Sets.difference(new HashSet<>(donorNames), refDonorNames);

		// validate that you have all the donor names in the reference data set.
		if (difference.size()>0) {
			log.error("Reference covariate file does not contain all donors in the covariate file.  Missing donors - " + difference.toString());
			return null;
		}
		
		// validate that you have all the requested covariates in the reference data set
		Set<String> refCovars=dc.getCovariates();
		SetView<String> differenceCovars = com.google.common.collect.Sets.difference(new HashSet<>(covars), refCovars);
		if (differenceCovars.size()>0) {
			log.error("Reference covariate file does not contain all requested covariates in the covariate file.  Missing covariates - " + differenceCovars.toString());
			return null;
		}

		for (String covar: covars) {
			String [] vals=dc.getValues(covar, donorNames).stream().toArray(String []::new);
			result.setValues(covar, vals);
		}
		return result;
	}
	
	private DonorCovariates parse(final List<File> covarFiles) {
		DonorCovariates result = DonorCovariates.parseFile(covarFiles.getFirst(), COVARIATE_VALIDATION);
		for (int i=1; i<covarFiles.size(); i++) {
			DonorCovariates other = DonorCovariates.parseFile(covarFiles.get(i), COVARIATE_VALIDATION);
			result.merge(other);
		}
		return (result);
	}
	
	

	public EqtlCovariate mergeCovariates(final List<File> covariateFiles, final File rejectedDonors, final boolean addBatchEffect) {
		Set<String> donors=Collections.emptySet();
		if (rejectedDonors!=null)
			donors = new HashSet<> (ParseBarcodeFile.readCellBarcodeFile(rejectedDonors));
		EqtlCovariate result = EqtlCovariate.parseFile(covariateFiles.getFirst(), donors, COVARIATE_VALIDATION);
		if (addBatchEffect) addBatchEffect(result, 1);
		for (int i=1; i<covariateFiles.size(); i++) {
			EqtlCovariate r = EqtlCovariate.parseFile(covariateFiles.get(i), donors, COVARIATE_VALIDATION);
			if (addBatchEffect) r = addBatchEffect(r, i+1);
			result.mergeDonors(r);
		}

		return result;
	}

	private EqtlCovariate addBatchEffect (final EqtlCovariate data, final int batchNum) {
		//String value = "batch"+batchNum;
		String value = Integer.toString(batchNum);
		String [] values = new String [data.donorNames().size()];
		Arrays.fill(values, value);
		data.setValues(BATCH_EFFECT_ATTRIBUTE, values);
		return data;
	}

	private EqtlCovariate addLibrarySize (final EqtlCovariate data, final File metaCellFile) {
		DGEMatrix m = DGEMatrix.parseDenseFile(metaCellFile, "");
		double [] totalExpression = m.getTotalExpressionPerSample();
		DecimalFormat format = new DecimalFormat("#");
		String [] values = Arrays.stream(totalExpression).mapToObj(x-> format.format(x)).toArray(String[]::new);
		data.setValues(this.LIBRARY_SIZE_ATTRIBUTE, values);
		return data;
	}

	public static EqtlCovariate generateEmptyCovariatesFile (final File metaCellFile) {
		// DGEMatrix m = DGEMatrix.parseDenseFile(f, prefix);
		String [] line = parseMetaCellHeader(metaCellFile, "GENE");
		List<String> donors = new ArrayList<> (Arrays.asList(line));
		// get rid of the first column, "GENE".
		donors.remove(0);
		EqtlCovariate result =new EqtlCovariate(donors);
		return result;

	}

	/**
	 * For parsing meta cell headers only.
	 * @param f
	 * @param expectedFirstColumn
	 * @return
	 */
	private static String [] parseMetaCellHeader (final File f, final String expectedFirstColumn) {
        TabbedInputParser parser = new TabbedInputParser(false, f);
        // read in lines until you're past any comment lines.
        String [] line=null;
        while (parser.hasNext()) {
        	line = parser.next();
        	if (line.length>0 && !line[0].startsWith("#")) break;
        }

        if (line==null)
			throw new IllegalArgumentException ("No data found in file " + f.getAbsolutePath());

        if (line.length==0 || !line[0].equals(expectedFirstColumn))
        	throw new IllegalArgumentException ("No "+expectedFirstColumn+" column as first entry in header.  Is this a meta cell file?" + StringUtils.join(line, "\t") + f.getAbsolutePath());
        return line;
	}


	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new PrepareEqtlCovariates().instanceMain(args));
	}


}
