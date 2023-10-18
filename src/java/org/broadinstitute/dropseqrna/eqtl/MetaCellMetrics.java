package org.broadinstitute.dropseqrna.eqtl;

import java.io.File;
import java.io.PrintStream;
import java.math.BigDecimal;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;

import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.IOUtil;
import picard.util.TabbedInputParser;

/**
 * A collection of metrics gathered over a collection of meta cells
 * @author nemesh
 *
 */
public class MetaCellMetrics {

	private final ObjectCounter<String> cellsPerDonor;
	private final Map<String, Double> umisPerDonor;
	private final LinkedHashSet<String> donors;
	private static final String [] header = {"DONOR", "NUM_CELLS", "NUM_UMIS", "AVERAGE_UMIS_PER_CELL"};
	
	/**
	 * Instantiate class with a list of donors that provide ordering for requested outputs.
	 */
	public MetaCellMetrics(List<String> donors) {
		this.donors = new LinkedHashSet<>(donors);
		this.cellsPerDonor=new ObjectCounter<>();
		this.umisPerDonor=new HashMap<>();
	}
	
	/**
	 * Add a single donor to the donor list.
	 * @param donor The donor to add
	 */
	public void addDonor (String donor) {
		this.donors.add(donor);
	}
	
	public boolean hasDonor (String donor) {
		return donors.contains(donor);
	}
	
	/**
	 * Set the number of cells for this donor.
	 * If the donor key already exists, then increment the total count by this amount.
	 * @param donor The donor to update
	 * @param cellCountForDonor The number of cells to attribute to this donor
	 */
	public void addCellsPerDonor (String donor, int cellCountForDonor) {
		if (!donors.contains(donor))
			throw new IllegalArgumentException("Donor not registered in metric collection");
		this.cellsPerDonor.incrementByCount(donor, cellCountForDonor);
	}
	
	/**
	 * Set the total number of UMIs for this donor.
	 * If the donor key already exists, then increment the total count by this amount.
	 * @param donor The donor to update
	 * @param umiCount The number of UMIs to attribute to this donor
	 */
	public void addUMIsPerDonor (String donor, double umiCount) {
		if (!donors.contains(donor))
			throw new IllegalArgumentException("Donor not registered in metric collection");
		Double val = umisPerDonor.get(donor);
		if (val!=null) 
			umiCount+=val;
		this.umisPerDonor.put(donor, umiCount);
		
		
	}
	
	/**
	 * Get a list of donors registered for the metrics.
	 * @return A list of donors in the order they were added.
	 */
	public List<String> getDonors() {
		return new ArrayList<>(this.donors);
	}
	
	/**
	 * Return the number of cells that belong to this donor
	 * @param donor The donor name
	 * @return The number of cells assigned to this donor
	 */
	public int getCellCount (String donor) {
		if (!donors.contains(donor))
			throw new IllegalArgumentException("Donor not registered in metric collection");
		return cellsPerDonor.getCountForKey(donor);
	}

	/**
	 * Return the number of UMIs that belong to this donor
	 * @param donor The donor name
	 * @return The number of UMIs assigned to this donor
	 */
	public double getUmiCount (String donor) {
		if (!donors.contains(donor))
			throw new IllegalArgumentException("Donor not registered in metric collection");
		if (this.umisPerDonor.containsKey(donor))
			return this.umisPerDonor.get(donor);
		return 0d;
	}
	
	/**
	 * Derive the mean number of UMIs per cell.
	 * @param donor The donor name
	 * @return The number of UMIs / the number of cells.
	 */
	public double getAverageUMIsPerCell (String donor) {
		return getUmiCount(donor) / (double) getCellCount(donor);
	}
	
	public void merge (MetaCellMetrics other) {
		for (String d: other.getDonors()) {
			if (!this.hasDonor(d)) this.addDonor(d);
			this.addCellsPerDonor(d, other.getCellCount(d));
			this.addUMIsPerDonor(d, other.getUmiCount(d));
		}
	}
	
	/**
	 * Parse a MetaCellMetric file and return a MetaCellMetrics object.
	 * @param f The file to parse
	 * @return A MetaCellMetrics object
	 */
	public static MetaCellMetrics readFile (File f) {
		// parse in the metric file.
		TabbedInputParser parser = new TabbedInputParser(false, f);
		Iterator<String [] > iter = parser.iterator();
		if (!iter.hasNext())
			throw new IllegalArgumentException ("Input meta cell metrics file is empty" + f.getAbsolutePath());
			
		String [] parsedHeader =iter.next();
		// map header to position in array.
		Map<String, Integer> position = new HashMap<>();
		for (int i=0; i<parsedHeader.length; i++) {
			position.put(parsedHeader[i], i);
		}
				
		MetaCellMetrics result = new MetaCellMetrics(new ArrayList<>());
		// parse the lines of the file.
		while (iter.hasNext()) {
			String [] line = iter.next();
			String donor = line [position.get("DONOR")];
			int cellCountForDonor = new BigDecimal(line[position.get("NUM_CELLS")]).intValue();
			int umiCount = new BigDecimal(line[position.get("NUM_UMIS")]).intValue();			
			result.addDonor(donor);
			result.addUMIsPerDonor(donor, umiCount);
			result.addCellsPerDonor(donor, cellCountForDonor);						
		}
		CloserUtil.close(parser);
		return result;		
	}
	
	/**
	 * Writes the object to an output file.
	 * @param outFile The file to write to
	 */
	public void writeMetrics (File outFile, boolean integerFormat) {
		DecimalFormat df = new DecimalFormat("#.##");
		PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(outFile));		
		out.println(StringUtils.join(header, "\t"));
		for (String donor: this.getDonors()) {			
			// format UMI count as double or string
			String umiCount=Double.toString(this.getUmiCount(donor));
			if (integerFormat) 
				umiCount=Integer.toString((int) this.getUmiCount(donor));						
			String [] line = {donor, Integer.toString(this.getCellCount(donor)), umiCount, df.format(this.getAverageUMIsPerCell(donor))};
			out.println(StringUtils.join(line, "\t"));
		}
		CloserUtil.close(out);
	}
	
	
	
	
	
}
