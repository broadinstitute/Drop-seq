/*
 * MIT License
 *
 * Copyright 2018 Broad Institute
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

import java.io.BufferedInputStream;
import java.io.File;
import java.io.PrintStream;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.lang.StringUtils;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.dropseqrna.barnyard.ParseBarcodeFile;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeader;
import org.broadinstitute.dropseqrna.barnyard.digitalexpression.DgeHeaderCodec;
import org.broadinstitute.dropseqrna.cmdline.DropSeq;
import org.broadinstitute.dropseqrna.utils.RetainRemoveList;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import picard.cmdline.CommandLineProgram;
import picard.cmdline.StandardOptionDefinitions;
import picard.util.TabbedInputParser;
import picard.util.TabbedTextFileWithHeaderParser;
import picard.util.TabbedTextFileWithHeaderParser.Row;

@CommandLineProgramProperties(
        summary = "Create a meta-cell file from a DGE file and a mapping file that details which cell barcodes belong to particular donors",
        oneLineSummary = "Summarize expression per donor",
        programGroup = DropSeq.class
)

public class CreateMetaCells extends CommandLineProgram {

	private final Log log = Log.getInstance(CreateMetaCells.class);

	@Argument(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The input DGE file.")
	public File INPUT;

	@Argument(doc = "A file that maps the cell barcode to the appropriate donor.  " +
			"This file has 2 columns with a header: cell, bestSample.  (Column labels may be overridden.)",
			mutex = {"SINGLE_METACELL_LABEL"})
	public File DONOR_MAP;
	
	@Argument(shortName="S", doc="The UMI summary strategy to use.")
    public DonorMergeStrategy STRATEGY = DonorMergeStrategy.Sum;

	@Argument(doc = "A file that contains a list of cell barcodes to exclude in the output donor summary.  Inclusions override exclusions if a cell barcode is listed in both files.", optional=true)
	public File EXCLUDE_CELL_BARCODES_FILE;

	@Argument(doc = "A file that contains a list of cell barcodes to include in the output donor summary.  Inclusions override exclusions if a cell barcode is listed in both files.", optional=true)
	public File INCLUDE_CELL_BARCODES_FILE;
	
	@Argument (doc="A file containing cluster assignments.  This assigns cell barcodes to groups.  File header columns: CELL_BARCODE, CLUSTER.  ", optional=true)
	public File CLUSTER_ASSIGNMENTS_FILE;
	
	@Argument (doc="A list of one of more cluster assignments to select cells from.  When provided, this finds the list of cell barcodes that have the same label in"
			+ "the CLUSTER_ASSIGNMENTS_FILE.  If CLUSTER_ASSIGNMENTS_FILE is provided, there must be at least one CLUSTER_ASSIGNMENT value provided.  "
			+ "The cell barcode labels may include a prefix string to disambiguate experiments that share the same cell barcode string, "
			+ "for example P60ENTSTNRep1P1_ACGT and P60ENTSTNRep2P1_ACGT.  The proper prefix can be supplied individually by the PREFIX argument, "
			+ "or can be supplied for a number of experiments using the MERGED_DGE_HEADER_FILE argument. ", optional=true)
	public List<String> CLUSTER_ASSIGNMENT;
	
	@Argument (doc="When running ICA Cluster, the merged DGE header contains both the UEI and PREFIX for each individual (DGE file) experiment.  "
			+ "This provides a way to look up the cell barcode prefix for this DGE file.  If using CLUSTER_ASSIGNMENTS_FILE, this or PREFIX must be set.", optional=true)
	public File MERGED_DGE_HEADER_FILE;
	
	@Argument (doc="If the CLUSTER_ASSIGNMENTS_FILE didn't come from ICA clustering, or there's some other workflow involved, you can provide the prefix for cell barcodes"
			+ "contained in the CLUSTER_ASSIGNMENTS_FILE directly.  If using CLUSTER_ASSIGNMENTS_FILE, this or MERGED_DGE_HEADER_FILE must be set.", optional=true)
	public String PREFIX;
	
	@Argument(shortName = StandardOptionDefinitions.OUTPUT_SHORT_NAME, doc="The output meta cell file.")
	public File OUTPUT;
	
	@Argument(shortName = StandardOptionDefinitions.METRICS_FILE_SHORT_NAME, doc="Metrics for each donor in the meta cell output.", optional=true)
	public File METRICS;
	
	@Argument(doc="Should the output DGE be formatted as integers?  Set this to true when merging normal DGE files.  Set to false "
			+ "if you're merging DGE files that you have performed some other transformations that change the expression results from "
			+ "quantized integers to a continuous variable")
	public Boolean INTEGER_FORMAT=true;

	@Argument(doc="Instead of creating multiple metacells, create a single metacell with this label.",
	mutex = {"DONOR_MAP", "CELL_BARCODE_COLUMN", "METACELL_COLUMN"})
	public String SINGLE_METACELL_LABEL;

	@Argument(doc="Column label for cell barcode column in DONOR_MAP.", mutex = {"SINGLE_METACELL_LABEL"})
	public String CELL_BARCODE_COLUMN = "cell";

	@Argument(doc="Column label for metacell column in DONOR_MAP.", mutex = {"SINGLE_METACELL_LABEL"})
	public String METACELL_COLUMN = "bestSample";

	private final String GENE_HEADER="GENE";
	private final int PROGRESS_INTERVAL=1000;
	
	private void validateClusterAssignmentArguments () {
		if (this.CLUSTER_ASSIGNMENTS_FILE!=null) {
			IOUtil.assertFileIsReadable(this.CLUSTER_ASSIGNMENTS_FILE);
			if (CLUSTER_ASSIGNMENT==null || CLUSTER_ASSIGNMENT.size()==0)
				throw new IllegalArgumentException("If CLUSTER_ASSIGNMENTS_FILE is set, then the CLUSTER_ASSIGNMENT must have at least one value.");
			if (MERGED_DGE_HEADER_FILE==null && PREFIX==null)
				throw new IllegalArgumentException("If CLUSTER_ASSIGNMENTS_FILE is set, then must set either MERGED_DGE_HEADER_FILE or PREFIX");
		}		
	}
	
	@Override
	public int doWork() {

		IOUtil.assertFileIsReadable(this.INPUT);
		if (this.DONOR_MAP != null) IOUtil.assertFileIsReadable(this.DONOR_MAP);
		IOUtil.assertFileIsWritable(this.OUTPUT);
		if (METRICS!=null) IOUtil.assertFileIsWritable(this.METRICS);
		if (INCLUDE_CELL_BARCODES_FILE!=null) IOUtil.assertFileIsReadable(this.INCLUDE_CELL_BARCODES_FILE);
		if (EXCLUDE_CELL_BARCODES_FILE!=null) IOUtil.assertFileIsReadable(this.EXCLUDE_CELL_BARCODES_FILE);
		validateClusterAssignmentArguments();
		
		PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(OUTPUT));

		final BufferedInputStream in = new BufferedInputStream(IOUtil.openFileForReading(INPUT));

		TabbedInputParser parser = new TabbedInputParser(false, in);
		if (!parser.hasNext()) {
			log.error("No lines in input file" + this.INPUT);
			parser.close();
			CloserUtil.close(out);
			return 1;
		}

		String [] header = parser.next();

		// what barcodes are in the output?  The first entry in the array is the gene name, so remove that.
		String [] cells = Arrays.copyOfRange(header, 1, header.length);
		// this is the list of cells to iterate over. calculate this just once.
		Set <String> cellsRetained = getElementsToRetain(Arrays.asList(cells), this.EXCLUDE_CELL_BARCODES_FILE, this.INCLUDE_CELL_BARCODES_FILE);
		log.info("Retaining [", cellsRetained.size(), "] cells");

		// read in how cells map to donors.
		Map<String,String> donorMap;
		final List<String> donors;
		if (this.DONOR_MAP != null) {
			donorMap = readDonorMap(this.DONOR_MAP, cellsRetained);
			
			// if there are clusters to subset data on, then subset the donor map by the cluster elements.
			
			String thisUEI = getUEI(this.INPUT);
			String prefix=getPrefixForExperiment (this.MERGED_DGE_HEADER_FILE, this.PREFIX, thisUEI);
			donorMap = filterDonorMapByClusterLabels (prefix, this.CLUSTER_ASSIGNMENT, this.CLUSTER_ASSIGNMENTS_FILE, donorMap);
			
			// the output order of the donors is the alpha numeric ordering of donor names.
			donors = new ArrayList<>(new HashSet<>(donorMap.values()));
			Collections.sort(donors);
		} else {
			donorMap = new HashMap<>();
			for (final String cell : cellsRetained) {
				donorMap.put(cell, SINGLE_METACELL_LABEL);
			}
			donors = Collections.singletonList(SINGLE_METACELL_LABEL);
		}
		
		// build meta cell statistics from donorMap
		MetaCellMetrics metrics = new MetaCellMetrics(donors);
		metrics = addCellCountsPerDonor(donorMap, metrics);
		
		// write the header for the output file.
		writeHeader(out, donors);

		// parse through the DGE and write out the meta-cell file.
		int counter=0;
		while (parser.hasNext()) {
			String [] line = parser.next();
			String gene = line[0];
			Map<String, Double> umis = getUMIsPerDonor (header, line, donorMap, this.STRATEGY);
			metrics = addUMIsPerDonor(metrics, umis);
			List<String> lineParsed = buildBodyLine(gene, umis, donors, this.INTEGER_FORMAT);
			if (lineParsed!=null)
				writeLine(lineParsed, out, false);
			if (counter > 0 && counter%PROGRESS_INTERVAL==0) log.info("Processed [" + counter +"] lines");
			counter++;
		}
		log.info("Processed [" + counter +"] lines in total");
		// cleanup.
		parser.close();
		out.close();
		
		// write metrics output if requested
		if (this.METRICS!=null) metrics.writeMetrics(this.METRICS, this.INTEGER_FORMAT); 			
		
		return 0;
	}

	String getPrefixForExperiment (final File mergedDgeHeader, final String userPrefix, final String uei) {
		if (userPrefix!=null) return userPrefix;
		if (mergedDgeHeader==null || uei==null) return (null);
		// if not declared up front, try the merged DGE header and look up matching UEI.
		DgeHeader mergedHeader = getDgeHeader(mergedDgeHeader);
		String resultPrefix=null;
		// iterate through the prefixes of the merged header, find the same UEI as the requested UEI. 
		Iterator<String> iter = mergedHeader.iterateLibraries();
		while (iter.hasNext()) {
			String mergedPrefix = iter.next();
			String mergedUEI = mergedHeader.getLibrary(mergedPrefix).getUei();
			if (mergedUEI.equals(uei))
				resultPrefix=mergedPrefix;
			}				
		return resultPrefix;
	}
	
	List<String> getCellBarcodesInClusters (String prefix, Set<String> clusterLabels, final File clusterAssignmentFile) {
		if (prefix==null || clusterAssignmentFile==null) return null;
		TabbedTextFileWithHeaderParser p = new TabbedTextFileWithHeaderParser(clusterAssignmentFile);
		Iterator<Row> iter =  p.iterator();
		
		List<String> result = new ArrayList<>();
		while (iter.hasNext()) {
			Row r = iter.next();
			String [] line = r.getFields();
			String cellBarcode = line[0];
			String cluster = line[1];
			int lastIndex = cellBarcode.lastIndexOf("_");
			String thisPrefix=cellBarcode.substring(0, lastIndex);
			String thisCellbarcode= cellBarcode.substring(lastIndex+1, cellBarcode.length());
			//String [] cellBCSplit = StringUtils.split(cellBarcode, "_");
			if (thisPrefix.equals(prefix) & clusterLabels.contains(cluster)) {
				result.add(thisCellbarcode);
			}			
		}
		return result;
	}
	
	private Map<String, String> filterDonorMapByClusterLabels (String prefix, Collection<String> clusterLabels, final File clusterAssignmentFile, final Map<String,String> donorMap) {
		// no op if prefix is unset or cluster assignment is null.
		if (prefix==null || clusterAssignmentFile==null) return (donorMap);
		if (clusterLabels==null || clusterLabels.size()==0) return (donorMap);
		
		Set<String> labels=new HashSet<>(clusterLabels);		
		Set<String> cb = new HashSet<> (getCellBarcodesInClusters(prefix, labels, clusterAssignmentFile));
		Map<String, String> result = new HashMap<>();
		for (String c: cb) {
			String v = donorMap.get(c);
			if (v!=null)
				result.put(c, v);
		}
		return result;
	}
	
	/**
	 * Get UEI for a DGE header that has a single library.
	 * @param input An input DGE file.
	 * @return The UEI if there is 1 and only 1 UEI in the file, else null.
	 */
	String getUEI (File input) {
		DgeHeader h = getDgeHeader(input);
		if (h==null) return null;
		int numLibs = h.getNumLibraries();
		if (numLibs==0) return null;
		if (numLibs!=1) 
			throw new IllegalArgumentException("Multiple UEIs found for file [" + input.getAbsolutePath() +"] are you sure this is a proper DGE file?");
		
		String uei = h.getLibrary(0).getUei();		
		return uei;		
	}
	
	
	private DgeHeader getDgeHeader (File input) {
		final BufferedInputStream inputStream = new BufferedInputStream(IOUtil.openFileForReading(input));
	    return new DgeHeaderCodec().decode(inputStream, input.getAbsolutePath());	    
	}
	
		
	MetaCellMetrics addCellCountsPerDonor (final Map<String,String> donorMap, MetaCellMetrics metrics) {
		donorMap.values().stream().forEach(x -> metrics.addCellsPerDonor(x, 1));
		return metrics;
	}
	
	MetaCellMetrics addUMIsPerDonor (MetaCellMetrics metrics, Map<String, Double> umis) {
		umis.keySet().forEach(x -> metrics.addUMIsPerDonor(x, umis.get(x)));
		return metrics;
	}

	/*
	private Map<String, Double> getUMIsPerDonorOld (String [] header, String [] line, final Map<String,String> donorMap) {
		// collapse the cell barcode data for this gene by donor name.
		// this should be a map to a list of doubles, which should then be summarized.
		Map<String, Double> umisPerDonor = new HashMap<>();
		for (int i=1; i<line.length; i++) {
			String cell = header[i];
			double count = Double.parseDouble(line[i]);
			String donor = donorMap.get(cell);
			if (donor!=null) {
				Double val = umisPerDonor.get(donor);
				if (val!=null) 
					count+=val;
				umisPerDonor.put(donor, count);					
			}				
		}
		return umisPerDonor;
	}
	*/
	
	private Map<String, Double> getUMIsPerDonor (String [] header, String [] line, final Map<String,String> donorMap, DonorMergeStrategy strategy) {
		// collapse the cell barcode data for this gene by donor name.
	
		// accumulate values per donor for this gene.
		Map<String, List<Double>> umiListPerDonor = new HashMap<>();
		
		for (int i=1; i<line.length; i++) {
			String cell = header[i];
			double count = Double.parseDouble(line[i]);
			String donor = donorMap.get(cell);			
			if (donor!=null) {
				List<Double> vals = umiListPerDonor.get(donor); 
				if (vals==null) {
					vals= new ArrayList<>();
					umiListPerDonor.put(donor, vals);
				}
				vals.add(count);					
			}				
		}		
		
		// finished accumulation, summarize by the selected strategy.		
		Map<String, Double> umisPerDonor = new HashMap<>();
		for (String donor: umiListPerDonor.keySet()) {
			List<Double> vals = umiListPerDonor.get(donor);
			double summarizedValue = summarizeDonor(vals, strategy);
			umisPerDonor.put(donor, summarizedValue);
		}
		
		return umisPerDonor;
	}
	
	private List<String> buildBodyLine (final String gene, Map<String,Double> umisPerDonor, List<String> donors, final boolean formatAsInteger) {
		// check if any values > 0.  If so, return the entire gene.  If not, return null.
		double totalCount = umisPerDonor.values().stream().mapToDouble(Double::doubleValue).sum();
		if (BigDecimal.ZERO.equals(new BigDecimal(totalCount))) return null;		
		
		// construct the output string.
		List<String> result = new ArrayList<>();		
		result.add(gene);
		for (String d: donors) {
			double v = umisPerDonor.get(d);			
			result.add(formatExpressionValue(v, formatAsInteger));
		}
		return result;
	}


	private Map<String,String> readDonorMap (final File donorMapFile, final Set <String> cellsRetained) {
		TabbedTextFileWithHeaderParser parser = new TabbedTextFileWithHeaderParser(donorMapFile);
		Map<String,String> result = new HashMap<>();
		for (final TabbedTextFileWithHeaderParser.Row row : parser) {
			final String cellBarcode = row.getField(CELL_BARCODE_COLUMN);
			if (cellsRetained.contains(cellBarcode)) {
				result.put(cellBarcode, row.getField(METACELL_COLUMN));
			}
		}
		CloserUtil.close (parser);
		return result;
	}


	/**
	 * Selects the set of elements(cells) to retain.
	 * @param elements The original list of cells
	 * @param remove A file containing a list of cells to remove, can be null.
	 * @param retain A file containing a list of cells to retain, can be null.
	 * @return The set of cell barcodes to retain.
	 */
	private Set<String> getElementsToRetain(final List<String> elements, final File remove, final File retain) {

		// get final set to remove
		Set <String> toRemove = getRetainRemoveSet(remove);
		Set <String> toRetain = getRetainRemoveSet(retain);

		// in the case where the input retain file is not null, but has 0 entries, there is nothing to retain. 
		if (retain!=null & toRetain.size()==0) 
			return Collections.emptySet();		

		RetainRemoveList<String> rrl = new RetainRemoveList<>();
		List<String> t =rrl.getElementsToRetain(elements, toRemove, toRetain);
		return new HashSet<> (t);

	}

	private Set<String> getRetainRemoveSet (final File f) {
		if (f==null)
			return (new HashSet<>(0));
		return (new HashSet<> (ParseBarcodeFile.readCellBarcodeFile(f)));
	}

	private void writeHeader (final PrintStream out, final List<String> donors) {
		out.println("#"+getCommandLine());
		List<String> lineOut = new ArrayList<>();

		// short-circuit: only write the gene header if there is at least one donor.
		if (donors.size()==0) return;
		
		lineOut.add(this.GENE_HEADER);
		lineOut.addAll(donors);
		String b = StringUtils.join(lineOut, "\t");
		out.println(b);
	}
	
	private void writeLine (final List<String> line, final PrintStream out, final boolean addGeneHeader) {
		List<String> lineOut = new ArrayList<>(line);
		if (addGeneHeader)
			lineOut.add(0, this.GENE_HEADER);
		String b = StringUtils.join(lineOut, "\t");
		out.println(b);
	}
	
	private String formatExpressionValue(final double exp, final boolean formatAsInteger) {
		if (formatAsInteger)
			return Integer.toString((int)exp);
		else
			return Double.toString(exp);
	}
	
	private Double summarizeDonor (final List<Double> expression, DonorMergeStrategy strategy) {
		if (strategy==DonorMergeStrategy.Sum) {
			return expression.stream().mapToDouble(Double::doubleValue).sum();
		}
		if (strategy==DonorMergeStrategy.Mean) {
			return expression.stream().mapToDouble(Double::doubleValue).average().getAsDouble();
		}
		if (strategy==DonorMergeStrategy.Median) {
			Median m = new Median();
			double [] values = expression.stream().mapToDouble(Double::doubleValue).toArray();
			return m.evaluate(values);
		}
		return null;
	}
	
	 


	/** Stock main method. */
	public static void main(final String[] args) {
		System.exit(new CreateMetaCells().instanceMain(args));
	}


}
