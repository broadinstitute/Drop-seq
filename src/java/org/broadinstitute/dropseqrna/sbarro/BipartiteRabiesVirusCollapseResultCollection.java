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

package org.broadinstitute.dropseqrna.sbarro;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;

/**
 * Structure to hold a single rabies virus, a list of potential intended sequences, the edit distances, and the UMI sharing.
 * Each intended sequence is a collapse result.
 * The collection only holds multiple results if there are multiple results where a tie needs to be broken (multiple parents within edit distance, etc.)
 * @author nemesh
 *
 */
public class BipartiteRabiesVirusCollapseResultCollection {

	private final List<BipartiteRabiesVirusCollapseResult> results;
	private BipartiteRabiesVirusCollapseResult bestResult=null;
	private final String childBarcode;
	private final static double EPSILON=0.0001;

	public BipartiteRabiesVirusCollapseResultCollection (final String childBarcode) {
		this.childBarcode=childBarcode;
		results=new ArrayList<>();
	}

	public void add (final String parentBarcode, final int editDistanceLeft, final int editDistanceRight, final double umiSharing) {
		BipartiteRabiesVirusCollapseResult r = new BipartiteRabiesVirusCollapseResult(this.childBarcode, parentBarcode, editDistanceLeft, editDistanceRight, umiSharing);
		this.results.add(r);
	}

	public void add (final String parentBarcode, final int editDistanceLeft, final int editDistanceRight) {
		BipartiteRabiesVirusCollapseResult r = new BipartiteRabiesVirusCollapseResult(this.childBarcode, parentBarcode, editDistanceLeft, editDistanceRight);
		this.results.add(r);
	}

	public void add (final BipartiteRabiesVirusCollapseResult r) {
		if (!r.getChildBarcode().equals(this.childBarcode))
			throw new IllegalArgumentException("Trying to add a collapse result with child sequence " + r.getChildBarcode() + " when the collection is for sequence " + this.childBarcode);
		this.results.add(r);
	}

	public List<BipartiteRabiesVirusCollapseResult> getResults() {
		return results;
	}

	public String getChildBarcode() {
		return childBarcode;
	}

	/**
	 * Get a map from each small barcode to the related larger barcodes.
	 * @return
	 */
	public Map<String,Collection<String>> getBarcodeMapping () {
		Multimap<String, String> result = ArrayListMultimap.create();
		for (BipartiteRabiesVirusCollapseResult b: this.results)
			result.put(b.getChildBarcode(), b.getParentBarcode());
		return result.asMap();
	}

	public BipartiteRabiesVirusCollapseResult getBestResult() {
		// short circuit for speed if this is called more than once.
		if (this.bestResult!=null) return this.bestResult;

		double bestUMISharing=Double.MIN_VALUE;
		List<BipartiteRabiesVirusCollapseResult> bestResults = new ArrayList<>();

		for (BipartiteRabiesVirusCollapseResult r: results) {
			Double umiSharing = r.getUmiSharing();
			if (umiSharing!=null) {
				if (Math.abs(umiSharing- bestUMISharing) < EPSILON)
					bestResults.add(r);
				else if (umiSharing>bestUMISharing) {
					bestUMISharing=umiSharing;
					bestResults.clear();
					bestResults.add(r);
				}
			} else // if there's no umi sharing data, everything passes.
				bestResults.add(r);
		}

		// break ties in umi sharing with parent size.
		// if you don't have anything here, the first result is the best.
		BipartiteRabiesVirusCollapseResult bestResult = bestResults.get(0);
		if (bestResults.size()>1) {
			int bestUmiCount=Integer.MIN_VALUE;
			for (BipartiteRabiesVirusCollapseResult r: bestResults) {
				Integer umiCount = r.getUmisParent();
				if (umiCount!=null && umiCount>bestUmiCount) {
					bestResult=r;
					bestUmiCount=umiCount;
				}
			}
		}
		this.bestResult=bestResult;
		return bestResult;
	}

	@Override
	public String toString () {
		if (this.bestResult==null) getBestResult();
		return this.bestResult.toString();
	}



}
