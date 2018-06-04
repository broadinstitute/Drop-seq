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

package org.broadinstitute.dropseqrna.beadsynthesis;

import org.broadinstitute.dropseqrna.utils.ObjectCounter;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.HashMap;
import java.util.Map;



public class IntendedSequenceBuilderTest {

	@Test
	public void testIndel () {
		String intendedSequence ="AAAAACCGGGAA";
		BarcodeNeighborGroup neighbors = new BarcodeNeighborGroup("AAAACCGGGAAN");
		neighbors.addNeighbor(new BeadSynthesisErrorData("AAAACCGGGAAA"));
		neighbors.addNeighbor(new BeadSynthesisErrorData("AAAACCGGGAAC"));
		neighbors.addNeighbor(new BeadSynthesisErrorData("AAAACCGGGAAG"));
		neighbors.addNeighbor(new BeadSynthesisErrorData("AAAACCGGGAAT"));

		ObjectCounter<String> umiCounts = new ObjectCounter<>();
		umiCounts.incrementByCount("AAAAACCGGGAA", 6130);
		umiCounts.incrementByCount("AAAACCGGGAAA", 60);
		umiCounts.incrementByCount("AAAACCGGGAAC", 60);
		umiCounts.incrementByCount("AAAACCGGGAAG", 60);
		umiCounts.incrementByCount("AAAACCGGGAAT", 60);

		Map<String, Double> umiBias = new HashMap<>();
		umiBias.put("AAAAACCGGGAA", 0.008);
		umiBias.put("AAAACCGGGAAA", 0.9845);
		umiBias.put("AAAACCGGGAAC", 0.9845);
		umiBias.put("AAAACCGGGAAG", 0.9845);
		umiBias.put("AAAACCGGGAAT", 0.9845);

		IntendedSequenceBuilder b = new IntendedSequenceBuilder(umiCounts, umiBias);
		IntendedSequence result = b.build(intendedSequence, neighbors);

		Assert.assertEquals(result.getIntendedSequence(), "AAAAACCGGGAA");
		Assert.assertEquals(result.getDeletedBase(), new Character ('A'));
		Assert.assertEquals(result.getDeletedBasePos(), new Integer (5));
		Assert.assertEquals(result.getDeletionRate(), 0.03767661, 0.001);
		Assert.assertEquals(result.getIntendedSequenceUMIBias(), 0.008, 0.001);
		Assert.assertEquals(result.getRelatedMedianUMIBias(), 0.9845, 0.001);
		Assert.assertEquals(result.getIntendedSequenceUMIs(), new Integer(6130));
		Assert.assertEquals(result.getMedianRelatedSequenceUMIs(), 60, 0.01);
	}
	/** group  intendedSeq                       relatedSequences deletedBase
 	TGTATTGTTGG TGTATTGTTGGG TGTATTGTTGGA:TGTATTGTTGGC:TGTATTGTTGGT        <NA>
    deletedBasePos rate intendedTBias medianNeighborTBias intendedUMIs
             	NA   NA         0.963                   1           54
    neighborMaxUMIs neighborMedianUMIs neighborUMICOV neighborMinTBias
              	51                 49     0.03117398             0.98
	*/
	// in this case the rate, basepos, base deleted are null.  There's no intended sequence, but there is a "root" sequence.

	@Test
	public void fullDeletion () {
		String intendedSequence =null;
		BarcodeNeighborGroup neighbors = new BarcodeNeighborGroup("TGTATTGTTGGN");
		neighbors.addNeighbor(new BeadSynthesisErrorData("TGTATTGTTGGA"));
		neighbors.addNeighbor(new BeadSynthesisErrorData("TGTATTGTTGGC"));
		neighbors.addNeighbor(new BeadSynthesisErrorData("TGTATTGTTGGG"));
		neighbors.addNeighbor(new BeadSynthesisErrorData("TGTATTGTTGGT"));

		ObjectCounter<String> umiCounts = new ObjectCounter<>();
		umiCounts.incrementByCount("TGTATTGTTGGA", 54);
		umiCounts.incrementByCount("TGTATTGTTGGC", 49);
		umiCounts.incrementByCount("TGTATTGTTGGG", 51);
		umiCounts.incrementByCount("TGTATTGTTGGT", 50);

		Map<String, Double> umiBias = new HashMap<>();
		umiBias.put("TGTATTGTTGGA", 0.98);
		umiBias.put("TGTATTGTTGGC", 0.98);
		umiBias.put("TGTATTGTTGGG", 0.963);
		umiBias.put("TGTATTGTTGGT", 0.98);


		IntendedSequenceBuilder b = new IntendedSequenceBuilder(umiCounts, umiBias);
		IntendedSequence result = b.build(intendedSequence, neighbors);

		Assert.assertNull(result.getIntendedSequence());
		Assert.assertNull(result.getDeletedBase());
		Assert.assertNull(result.getDeletedBasePos());
		Assert.assertNull(result.getDeletionRate());
		Assert.assertNull(result.getIntendedSequenceUMIBias());
		Assert.assertEquals(result.getRelatedMedianUMIBias(), 0.98, 0.001);
		Assert.assertNull(result.getIntendedSequenceUMIs());
		Assert.assertEquals(result.getMedianRelatedSequenceUMIs(),50.5, 0.01);
	}

	@Test
	public void partialBase12Deletion () {
		String intendedSequence ="TGTATTGTTGGA";
		BarcodeNeighborGroup neighbors = new BarcodeNeighborGroup("TGTATTGTTGGN");
		neighbors.addNeighbor(new BeadSynthesisErrorData("TGTATTGTTGGC"));
		neighbors.addNeighbor(new BeadSynthesisErrorData("TGTATTGTTGGG"));
		neighbors.addNeighbor(new BeadSynthesisErrorData("TGTATTGTTGGT"));

		ObjectCounter<String> umiCounts = new ObjectCounter<>();
		umiCounts.incrementByCount("TGTATTGTTGGA", 50);
		umiCounts.incrementByCount("TGTATTGTTGGC", 50);
		umiCounts.incrementByCount("TGTATTGTTGGG", 50);
		umiCounts.incrementByCount("TGTATTGTTGGT", 50);

		Map<String, Double> umiBias = new HashMap<>();
		umiBias.put("TGTATTGTTGGA", 0.5);
		umiBias.put("TGTATTGTTGGC", 0.98);
		umiBias.put("TGTATTGTTGGG", 0.963);
		umiBias.put("TGTATTGTTGGT", 0.98);


		IntendedSequenceBuilder b = new IntendedSequenceBuilder(umiCounts, umiBias);
		IntendedSequence result = b.build(intendedSequence, neighbors);

		Assert.assertEquals(result.getIntendedSequence(), "TGTATTGTTGGA");
		Assert.assertEquals(result.getDeletedBase(), new Character ('A'));
		Assert.assertEquals(result.getDeletedBasePos(), new Integer (12));
		Assert.assertEquals(result.getDeletionRate(), 0.750, 0.001);
		Assert.assertEquals(result.getIntendedSequenceUMIBias(), 0.5, 0.001);
		Assert.assertEquals(result.getRelatedMedianUMIBias(), 0.98, 0.001);
		Assert.assertEquals(result.getIntendedSequenceUMIs(), new Integer(50));
		Assert.assertEquals(result.getMedianRelatedSequenceUMIs(), 50, 0.01);

	}


}
