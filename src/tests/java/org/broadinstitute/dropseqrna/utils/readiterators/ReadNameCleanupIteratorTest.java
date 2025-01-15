package org.broadinstitute.dropseqrna.utils.readiterators;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordSetBuilder;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Pattern;

import static org.testng.Assert.*;

public class ReadNameCleanupIteratorTest {

    private List<SAMRecord> records;
    private Iterator<SAMRecord> underlyingIterator;

    @BeforeMethod
    public void setUp() {
        SAMRecordSetBuilder builder = new SAMRecordSetBuilder();
        builder.addUnmappedFragment("prefix_read1_suffix");
        builder.addUnmappedFragment("prefix_read2_suffix");
        records = new ArrayList<>();
        builder.iterator().forEachRemaining(records::add);
    }

    @Test
    public void testPrefixRemoval() {
        underlyingIterator = records.iterator();
        ReadNameCleanupIterator iterator = new ReadNameCleanupIterator.Builder(underlyingIterator)
                .prefixToRemove("prefix_")
                .build();

        SAMRecord result = iterator.next();
        assertEquals(result.getReadName(), "read1_suffix");
    }

    @Test
    public void testSuffixRemoval() {
        underlyingIterator = records.iterator();
        ReadNameCleanupIterator iterator = new ReadNameCleanupIterator.Builder(underlyingIterator)
                .suffixToRemove("_suffix")
                .build();

        SAMRecord result = iterator.next();
        assertEquals(result.getReadName(), "prefix_read1");
    }

    @Test
    public void testPatternRemoval() {
        underlyingIterator = records.iterator();
        ReadNameCleanupIterator iterator = new ReadNameCleanupIterator.Builder(underlyingIterator)
                .patternToRemove(Pattern.compile("_read1_"))
                .build();

        SAMRecord result = iterator.next();
        assertEquals(result.getReadName(), "prefixsuffix");
    }

    @Test
    public void testPrefixAndSuffixAddition() {
        underlyingIterator = records.iterator();
        ReadNameCleanupIterator iterator = new ReadNameCleanupIterator.Builder(underlyingIterator)
                .prefixToAdd("new_prefix_")
                .suffixToAdd("_new_suffix")
                .build();

        SAMRecord result = iterator.next();
        assertEquals(result.getReadName(), "new_prefix_prefix_read1_suffix_new_suffix");
    }

    @Test
    public void testCombinedTransformations() {
        underlyingIterator = records.iterator();
        ReadNameCleanupIterator iterator = new ReadNameCleanupIterator.Builder(underlyingIterator)
                .prefixToRemove("prefix_")
                .suffixToRemove("_suffix")
                .prefixToAdd("new_prefix_")
                .suffixToAdd("_new_suffix")
                .build();

        SAMRecord result = iterator.next();
        assertEquals(result.getReadName(), "new_prefix_read1_new_suffix");
    }
}