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
package org.broadinstitute.dropseqrna.utils.modularfileparser;

import java.io.File;

import org.testng.annotations.Test;

import org.testng.Assert;

public class ParserTest {

	private File testFile1= new File ("testdata/org/broadinstitute/transcriptome/utils/modularfileparser/ClozUK_CNV_Loci.txt");
    private File bedFile = new File ("testdata/org/broadinstitute/transcriptome/utils/modularfileparser/testBed.bed.txt");

	@Test()
	//TODO: these tests kinda suck.
	public void DelimiterParserTest() {

		Parser p = new DelimiterParser("\t");
		ModularFileParser mfp = null;
		mfp = new ModularFileParser(p, testFile1, 0);
		String [] body = null;
		while ((body=mfp.readNextLine())!=null)
			Assert.assertNotNull(body);
	}

	@Test()
	//TODO: these tests kinda suck.
	public void BEDParserTest() {
		BEDFileParser p = new BEDFileParser(" ");
		ModularFileParser mfp = new ModularFileParser(p, bedFile, 1);
		String [] body = null;
		while ((body=mfp.readNextLine())!=null)
			Assert.assertNotNull(body);
		String [] header =p.getDefaultHeader();
		Assert.assertNotNull(header);
	}
}
