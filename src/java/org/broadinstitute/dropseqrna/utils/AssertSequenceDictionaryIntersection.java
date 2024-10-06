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
package org.broadinstitute.dropseqrna.utils;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.CloserUtil;
import htsjdk.samtools.util.Log;
import htsjdk.variant.vcf.VCFFileReader;

import java.io.File;
import java.nio.file.Path;

/**
 * Throw an exception if there is no overlap of sequence names in the two inputs.
 * If log is non-null, log a INFO message if there is an overlap.
 */
public class AssertSequenceDictionaryIntersection {

    public static void assertIntersectionVcfBam(final File vcf, final File bam, final Log log) {
        assertIntersectionVcfBam(vcf.toPath(), bam.toPath(), log);
    }

    public static void assertIntersectionVcfBam(final Path vcf, final Path bam, final Log log) {
        final VCFFileReader vcfReader = new VCFFileReader(vcf, false);
        try {
            assertIntersectionObjectBam(vcfReader, vcf.getFileName().toString(), bam, log);
        } finally {
            CloserUtil.close(vcfReader);
        }
    }

    /**
     * @param obj Some object from which a sequence dictionary can be extracted
     * @param objDescription user-friendly description of obj, or null
     */
    public static void assertIntersectionObjectBam(final Object obj, final String objDescription, final File bam, final Log log) {
        assertIntersectionObjectBam(obj, objDescription, bam.toPath(), log);
    }

    public static void assertIntersectionObjectBam(final Object obj, final String objDescription, final Path bam,
                                                   final Log log) {
        final SamReader samReader = SamReaderFactory.makeDefault().open(bam);
        try {
            assertIntersection(obj, objDescription, samReader, bam.getFileName().toString(), log);
        } finally {
            CloserUtil.close(samReader);
        }
    }

    /**
     * @param obj Some object from which a sequence dictionary can be extracted
     * @param objDescription user-friendly description of obj, or null
     */
    public static void assertIntersectionObjectVcf(final Object obj, final String objDescription, final File vcf, final Log log) {
        assertIntersectionObjectVcf(obj, objDescription, vcf.toPath(), log);
    }

    /**
     * @param obj Some object from which a sequence dictionary can be extracted
     * @param objDescription user-friendly description of obj, or null
     */
    public static void assertIntersectionObjectVcf(final Object obj, final String objDescription, final Path vcf,
                                                   final Log log) {
        final VCFFileReader vcfReader = new VCFFileReader(vcf, false);
        try {
            assertIntersection(obj, objDescription, vcfReader, vcf.getFileName().toString(), log);
        } finally {
            CloserUtil.close(vcfReader);
        }
    }

    /**
     * Both descriptions can be null, and will be replaced with a generic string based on data type.
     */
    public static void assertIntersection(final Object obj1, final String description1,
                                          final Object obj2, final String description2, final Log log) {
        final SequenceDictionaryIntersection sdi = new SequenceDictionaryIntersection(
                obj1, description1, obj2, description2);
        final String message = sdi.message(false);
        if (sdi.getIntersection().isEmpty()) {
            throw new RuntimeException(message);
        } else {
            if (log != null) {
                log.info(message);
            }
        }
    }
}
