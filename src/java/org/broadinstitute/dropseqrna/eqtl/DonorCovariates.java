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

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Sets.SetView;

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import org.broadinstitute.dropseqrna.eqtl.NonNumericCovariate;
import picard.util.TabbedInputParser;

/**
 * Holds the covariate values for some donors, parsed in from a file.
 * Each donor should have some value for each covariate ie : a dense matrix of covariates.
 * @author nemesh
 *
 */
public class DonorCovariates {

	private static final Log LOG = Log.getInstance(DonorCovariates.class);
	
	// key is the donor name
	// value is a map of covariates to their values for a donor.
	private Map<String, Map<String,String>> values;

	private DonorCovariates () {
		values=new HashMap<>();
	}

	public void addValue (final String donor, final String covariate, final String value) {
		Map<String,String> v = this.values.get(donor);
		if (v==null) {
			v=new HashMap<>();
			this.values.put(donor, v);
		}
		v.put(covariate, value);
	}
	
	/**
	 * Get the set of covariates defined for the donors
	 */
	public Set<String> getCovariates () {
		Set<String> result = new HashSet<String>();
		values.values().forEach(x -> result.addAll((Collection<? extends String>) x.keySet()));
		return result;
	}

	/**
	 * Gets the value of a covariate for a donor.
	 * @param donor
	 * @param covariate
	 * @return
	 */
	public String getValue(final String donor, final String covariate) {
		Map<String,String> v = this.values.get(donor);
		if (v==null) return null;
		String vv = v.get(covariate);
		return vv;
	}

	public List<String> getValues(final String covariate, final List<String>donors) {
		List<String> result = new ArrayList<>();
		for (String d: donors)
			result.add(getValue(d, covariate));
		return result;
	}

	public List<String> getDonorNames () {
		return new ArrayList<>(this.values.keySet());
	}
	
	/**
	 * Merges covariates from another data set into this one.
	 * Merging covariates where donors don't overlap in some sets will result in covariates that are missing for the setdiff set of donors.
	 * @param other Another set of covariates
	 */
	public void merge (DonorCovariates other) {
		SetView<String>  differentDonors = com.google.common.collect.Sets.difference(new HashSet<>(other.getDonorNames()), new HashSet<>(this.getDonorNames()));
		if (!differentDonors.isEmpty()) {
			LOG.warn(
					"When merging covariates from multiple sources, new donor names detected that may result in some missing values " + differentDonors
			);
		}
		for (String donor: other.getDonorNames()) {
			for (String covariate: other.getCovariates()) {
				String value = other.getValue(donor, covariate);
				this.addValue(donor, covariate, value);
			}
		}
 	}

	/**
	 * Parse a covariate input file.
	 * The file format is:
	 * Header with donor ID column followed by covariate names
	 * Each line has a donor ID followed by covariate values
	 * File is tab sep.
	 * @param input The file to parse.
	 * @param validationStringency How to handle non-numeric covariates.
	 * @return A DonorCovariates object.
	 */
	public static DonorCovariates parseFile(final File input, final ValidationStringency validationStringency) {
		final TabbedInputParser parser = new TabbedInputParser(false, input);
		DonorCovariates result = new DonorCovariates();
        if (!parser.hasNext()) {
        	CloserUtil.close(parser);
        	return result;
        }

        // parse header.
        String [] header = parser.next();
        // assert this is a valid covariate reference or metrics file.
        if (!(header[0].equals("IID") || header[0].equals("DONOR"))) {
        	CloserUtil.close(parser);
        	throw new IllegalArgumentException("Expected header first row/column should be either IID or DONOR");        	
        }
        int numColumns=header.length;
		final List<NonNumericCovariate> nonNumericCovariates = new ArrayList<>();

        while (parser.hasNext()){
        	String [] line = parser.next();
        	String donorName=line[0];

        	if (line.length != numColumns) {
        		CloserUtil.close(parser);
				throw new RuntimeException(String.format("On line %d, expected %d columns but saw %d",
                        parser.getCurrentLineNumber(), numColumns, line.length));
        	}
        	for (int i=1; i<line.length; i++) {
        		final String field = line[i];
				final String attribute = header[i];
        		result.addValue(donorName, attribute, field);
				try {
					Double.parseDouble(field);
				} catch (NumberFormatException e) {
					nonNumericCovariates.add(new NonNumericCovariate(input, attribute, donorName, field));
				}
			}
		}
		CloserUtil.close(parser);
		if (!nonNumericCovariates.isEmpty()) {
			final String message = NonNumericCovariate.generateNonNumericMessage(nonNumericCovariates);
			switch(validationStringency) {
				case STRICT:
					throw new IllegalArgumentException(message);
				case LENIENT:
					LOG.warn(message);
					break;
				case SILENT:
					break;
			}
		}
		return result;
	}
}
