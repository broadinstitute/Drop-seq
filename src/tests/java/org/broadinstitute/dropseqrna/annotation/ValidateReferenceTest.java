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
package org.broadinstitute.dropseqrna.annotation;

import com.google.gson.Gson;
import htsjdk.samtools.util.CloserUtil;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class ValidateReferenceTest {
    private static final File TEST_DATA_DIR = new File("testdata/org/broadinstitute/transcriptome/annotation/");

    @Test
    public void testTrivial() {
        final String refName = "ERCC92";
        final String[] args = new String[]{
                "GTF=" + new File(TEST_DATA_DIR, refName + ".gtf.gz").getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + new File(TEST_DATA_DIR, refName + ".fasta.gz").getAbsolutePath()
        };
        Assert.assertEquals(new ValidateReference().instanceMain(args), 0);
    }

    @Test
    public void testProblems() throws IOException {
        final String refName = "buggy";
        final File output = File.createTempFile("ValidateReferenceTest.", ".json");
        output.deleteOnExit();
        final String[] args = new String[]{
                "GTF=" + new File(TEST_DATA_DIR, refName + ".gtf").getAbsolutePath(),
                "REFERENCE_SEQUENCE=" + new File(TEST_DATA_DIR, refName + ".fasta").getAbsolutePath(),
                "OUTPUT=" + output.getAbsolutePath()
        };
        Assert.assertEquals(new ValidateReference().instanceMain(args), 0);
        final FileReader reader = new FileReader(output);
        final ValidateReference.Messages messages = new Gson().fromJson(reader, ValidateReference.Messages.class);
        CloserUtil.close(reader);
        Assert.assertEquals(messages.transcriptsWithNoExons.size(), 0);
        Assert.assertNotNull(messages.baseErrors);
        Assert.assertTrue(messages.baseErrors.contains("88"));
        Assert.assertEquals(messages.sequencesOnlyInReference.size(), 2);
        Assert.assertTrue(messages.sequencesOnlyInReference.contains("ERCC_00171"));
        Assert.assertTrue(messages.sequencesOnlyInReference.contains("NOT_IN_GTF"));
        Assert.assertEquals(messages.sequencesOnlyInGtf.size(), 1);
        Assert.assertEquals(messages.sequencesOnlyInGtf.get(0), "BOGUS_SEQUENCE");
        Assert.assertEquals(messages.geneBiotypes.size(), 1);
        Assert.assertTrue(messages.geneBiotypes.contains(null));
        Assert.assertEquals(messages.fractionOfSequencesOnlyInReference, 2/93.0, 0.001);
        Assert.assertEquals(messages.fractionOfSequencesOnlyInGtf, 1/92.0, 0.001);
        Assert.assertEquals(messages.fractionOfGenomeOfSequencesOnlyInReference, 551.0/82802, 0.0001);
    }
}
