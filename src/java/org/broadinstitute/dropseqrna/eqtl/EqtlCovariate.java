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

import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.PeekableIterator;
import htsjdk.samtools.util.StringUtil;
import org.broadinstitute.dropseqrna.TranscriptomeException;
import org.broadinstitute.dropseqrna.utils.io.ErrorCheckingPrintStream;

import com.google.common.collect.Sets;

import picard.util.TabbedInputParser;

import java.io.File;
import java.io.PrintStream;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * A eQTL covariate is a matrix of donor IDs in columns and attributes in rows.
 * A valid eQTL should have a value for each donor/attribute pairing - ie a complete matrix.
 * @author nemesh
 *
 */
public class EqtlCovariate {

	// attribute to list of donors values for that attribute.
	private Map<String, String []> data;
	private List<String> donorNames;

	// the 0th row/column has this string in every eQTL matrix.
	public static final String FIRST_ELEMENT="id";
	private static final String MISSING_VALUE = "NA";

	private static final Log LOG = Log.getInstance(EqtlCovariate.class);

	public EqtlCovariate (final Collection<String> donors) {
		data = new LinkedHashMap<>();
		donorNames = new ArrayList<>(donors);
	}

	
	public void setValues (final String attribute, final String [] values) {
		if (values.length!=this.donorNames().size())
			throw new IllegalArgumentException("Attempting to set values on attribute ["+ attribute +"] for the wrong number of donors.  [" + this.donorNames().size()+" donors in collection [" + values.length +"] values provided");
		// replace nulls with missing values
		for (int i=0; i<values.length; i++) {
			if (values[i]==null)
				values[i] = MISSING_VALUE;
		}
        data.putIfAbsent(attribute, values);
	}
		
	public String [] getValues (final String attribute) {
		String [] d= data.get(attribute);
		return d;
	}

	public Set<String> getAttributes () {
		return this.data.keySet();
	}

	public List<String> donorNames () {
		return this.donorNames;
	}

	/**
	 * EQTL objects must have the same sets of attributes to be valid for merging donors
	 * Merges the data from <OTHER> into this EqtlCovariate object.
	 * @param other
	 * @return
	 */
	public void mergeDonors (final EqtlCovariate other) {
		if (!validateMergeDonors(other))
			throw new TranscriptomeException("The eQTL covariate attributes don't match, they can not be merged. "
					+ " Attribute set 1 "+ this.getAttributes() + " attribute set 2 " + other.getAttributes());

		Collection<String> donors = other.donorNames();
		this.donorNames.addAll(donors);
		for (String attribute: other.getAttributes()) {
			// lambdas, still magic.
			String [] result = Stream.of(getValues(attribute), other.getValues(attribute)).flatMap(Stream::of).toArray(String[]::new);
			this.data.put(attribute, result);
		}
	}
	
	public void mergeAttributes (final EqtlCovariate other) {
		if (!validateMergeAttributes(other))
			throw new TranscriptomeException("The eQTL covariate donors don't match, they can not be merged. "
					+ " Donors set 1 "+ this.donorNames() + " donors set 2 " + other.donorNames());

		Set<String> attributesThis = getAttributes();
		Set<String> attributesOther = other.getAttributes();
		Set<String> newAttributes = Sets.difference(attributesOther, attributesThis);
		
		for (String attribute: newAttributes) {
			// lambdas, still magic.
			String [] result = other.getValues(attribute);
			this.data.put(attribute, result);
		}
	}
	

	/**
	 * The two eQTL objects must have the same set of attributes to be valid for merging donors.
	 * @param other
	 * @return
	 */
	public boolean validateMergeDonors (final EqtlCovariate other) {
		Set<String> attributesThis = getAttributes();
		Set<String> attributesOther = other.getAttributes();
		if (attributesThis.size()!=attributesOther.size()) return false;
		return (attributesThis.containsAll(attributesOther));
	}
	
	
	public boolean validateMergeAttributes (final EqtlCovariate other) {
		List<String> donorNamesThis = donorNames();
		List<String> donorNamesOther = other.donorNames();
		if (donorNamesThis.size()!=donorNamesOther.size()) return false;
		
		for (int i=0; i<donorNamesThis.size(); i++) {
			if (!donorNamesThis.get(i).equals(donorNamesOther.get(i)))
				return false;
		}
		return true;
	}

	public static EqtlCovariate parseFile(final File f,
										  final Set<String> ignoredDonors,
										  final ValidationStringency validationStringency) {
		IOUtil.assertFileIsReadable(f);
		@SuppressWarnings("resource")
		TabbedInputParser parser = new TabbedInputParser(false, f);
		PeekableIterator<String [] > iter = new PeekableIterator<>(parser.iterator());
		String  [] header = iter.next();
		if (header.length==0)
			throw new TranscriptomeException ("First non-comment line of eQTL Covariates matrix malformed - line empty.");
		if (!header[0].equals(EqtlCovariate.FIRST_ELEMENT))
			throw new TranscriptomeException ("Expected first line first column to be " + EqtlCovariate.FIRST_ELEMENT+ " found "+ header[0]);

		List<String> donors = new ArrayList<> (Arrays.asList(Arrays.copyOfRange(header, 1, header.length)));

		// find the columns to ignore.
		Set<Integer> colsToIgnore = new HashSet<>();
		for (int i=0; i<donors.size(); i++)
			if (ignoredDonors.contains(donors.get(i)))
				colsToIgnore.add(i);

		List<String> filteredDonors = new ArrayList<> (donors);
		filteredDonors.removeAll(ignoredDonors);

		final EqtlCovariate result = new EqtlCovariate(filteredDonors);
		final List<NonNumericCovariate> nonNumericCovariates = new ArrayList<>();

		while (iter.hasNext()) {
			String [] body = iter.next();
			if (body.length!=header.length)
				throw new TranscriptomeException(
						"Different number of columns in header and body :" +
								Arrays.toString(header) + " " + Arrays.toString(body)
				);
			final String attribute = body[0];

			// this is so awkward.  You want to only retain certain indexes in the array, but this is hard(er) to express cleanly in java compared to R/python array slices.
			String [] finalValues = new String [filteredDonors.size()];
			int idx=0;
			for (int i=1; i<body.length; i++)
				if (!colsToIgnore.contains(i-1)) {
					final String covariateValue = body[i];
					finalValues[idx]=covariateValue;
					try {
						Double.parseDouble(finalValues[idx]);
					} catch (final NullPointerException | NumberFormatException e) {
						nonNumericCovariates.add(
								new NonNumericCovariate(f, attribute, filteredDonors.get(idx), covariateValue)
						);
					}
                    idx++;
				}
			result.setValues(attribute, finalValues);
		}
		parser.close();

		if (!nonNumericCovariates.isEmpty()) {
			// make an error message with each covariate on a separate line
			final String message = NonNumericCovariate.generateNonNumericMessage(nonNumericCovariates);
			switch (validationStringency) {
				case STRICT:
					throw new TranscriptomeException(message);
				case LENIENT:
					LOG.warn(message);
					break;
				case SILENT:
					break;
			}
		}
		return result;
	}

	public static EqtlCovariate parseFile (final File f) {
		return parseFile(f, Collections.emptySet(), ValidationStringency.SILENT);
	}

	public void writeFile (final File f) {
		PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(f));
		List<String> h = new ArrayList<>();
		h.add(EqtlCovariate.FIRST_ELEMENT);
		h.addAll(this.donorNames());
		out.println(StringUtil.join("\t", h));
		for (String a: this.getAttributes()) {			
			List<String> body = new ArrayList<>();
			body.add(a);
			body.addAll(Arrays.asList(this.getValues(a)));
			out.println(StringUtil.join("\t", body));
		}
		out.close();
	}
	
	/**
	 * Write out the object to a file, with a desired ordering of donors.
	 * 
	 * @param f
	 * @param desiredDonorOrder
	 */
	public void writeFile (final File f, List<String> desiredDonorOrder) {
		// desiredDonorOrder can be a subset of all donors, or the donors in a different order
		List<String> finalDonorList=getOrderedIntersectDonors(desiredDonorOrder);
		// map the existing donors to their position in the data arrays
		Map<String, Integer> posMap = getDonorNameToPosition ();
		
		//find the ordering of donors in this object for the output.
		List<Integer> outPositions = new ArrayList<Integer>();
		for (String d: finalDonorList) {
			outPositions.add(posMap.get(d));
		}
				
		PrintStream out = new ErrorCheckingPrintStream(IOUtil.openFileForWriting(f));
		List<String> h = new ArrayList<>();
		h.add(EqtlCovariate.FIRST_ELEMENT);
		h.addAll(finalDonorList);
		out.println(StringUtil.join("\t", h));
		for (String a: this.getAttributes()) {			
			List<String> body = new ArrayList<>();
			body.add(a);
			String [] values = this.getValues(a);
			List<String> orderedValues = outPositions.stream().map(x-> values[x]).collect(Collectors.toList());
			body.addAll(orderedValues);
			out.println(StringUtil.join("\t", body));
		}
		out.close();
	}
	
	private List<String> getOrderedIntersectDonors (List<String> desiredDonorOrder) {
		List<String> result = new ArrayList<String> ();
		for (String d: desiredDonorOrder) {
			if (this.donorNames.contains(d))
				result.add(d);
		}
		return result;
	}
	
	private Map<String, Integer> getDonorNameToPosition () {
		Map<String, Integer> result = new HashMap<>();
		int counter=0;
		for (String d: this.donorNames) {
			result.put(d, counter);
			counter++;
		}
		return (result);
	}

	@Override
	public String toString () {
		StringBuilder b= new StringBuilder();
		// add header.
		List<String> h = new ArrayList<>();
		h.add(EqtlCovariate.FIRST_ELEMENT);
		h.addAll(this.donorNames());
		b.append(StringUtil.join("\t", h)+"\n");
		// add attributes.
		for (String a: this.getAttributes()) {
			List<String> body = new ArrayList<>();
			body.add(a);
			body.addAll(Arrays.asList(this.getValues(a)));
			b.append(StringUtil.join("\t", body)+"\n");
		}
		return b.toString();
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((data == null) ? 0 : data.hashCode());
		return result;
	}

	@Override
	public boolean equals(final Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		EqtlCovariate other = (EqtlCovariate) obj;
		if (data == null) {
			if (other.data != null)
				return false;
		} else if (!data.equals(other.data))
			return false;
		return true;
	}

}
