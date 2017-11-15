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
package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTagUtil;

import java.util.Comparator;

/**
 * Comparator for SAMRecord that orders by the value for the given tag, which must be a string value.
 * Not-present value for the tag sorts before any non-null value.
 */
public class StringTagComparator implements Comparator<SAMRecord> {
    private final short tag;

    public StringTagComparator(String tag) {
        this.tag = SAMTagUtil.getSingleton().makeBinaryTag(tag);
    }

    @Override
    public int compare(SAMRecord rec1, SAMRecord rec2) {
        final String s1 = (String)rec1.getAttribute(tag);
        final String s2 = (String)rec2.getAttribute(tag);

        if (s1 != null) {
            if (s2 == null)
                return 1;
            else {
                return s1.compareTo(s2);
            }
        } else if (s2 != null) {
            return -1;
        } else {
            return 0;
        }
    }
}
