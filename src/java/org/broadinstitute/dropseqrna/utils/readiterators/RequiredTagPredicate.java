/*
 * MIT License
 *
 * Copyright 2019 Broad Institute
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
package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTagUtil;

import java.util.function.Predicate;

public class RequiredTagPredicate
        implements Predicate<SAMRecord> {

    final short[] requiredTags;

    public RequiredTagPredicate(final String... requiredTags) {
        this.requiredTags = new short[requiredTags.length];
        for (int i = 0; i < requiredTags.length; ++i)
            this.requiredTags[i] = SAMTagUtil.getSingleton().makeBinaryTag(requiredTags[i]);
    }

    @Override
    public boolean test(SAMRecord rec) {
        for (final short tag : requiredTags)
            // String strTag=SAMTagUtil.getSingleton().makeStringTag(tag);
            // List<SAMTagAndValue> vals = rec.getAttributes();
            if (rec.getAttribute(tag) == null)
                return false;
        return true;
    }
}
