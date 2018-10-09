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

package org.broadinstitute.dropseqrna.utils;

import java.io.File;

import htsjdk.samtools.util.CloserUtil;
import picard.util.TabbedInputParser;

/**
 * Helper methods for
 * @author nemesh
 *
 */
public class TestUtils {

	/**
	 * Test if two text files are the same, ignoring "#" character lines.
	 * @param expected
	 * @param actual
	 * @return
	 */
	public static boolean testFilesSame (final File expected, final File actual) {
		TabbedInputParser e = new TabbedInputParser(true, expected);
		TabbedInputParser a = new TabbedInputParser(true, expected);
		while (e.hasNext() && a.hasNext()) {
			e.next();
			a.next();
			String le = e.getCurrentLine();
			String la = a.getCurrentLine();
			if (!le.equals(la)) {
				CloserUtil.close(e);
				CloserUtil.close(a);
				return false;
			}
		}
		CloserUtil.close(e);
		CloserUtil.close(a);
		return true;
	}
}
