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

public class BAMTagCleanupIteratorTest {

    private List<SAMRecord> records;
    private Iterator<SAMRecord> underlyingIterator;

    @BeforeMethod
    public void setUp() {
        SAMRecordSetBuilder builder = new SAMRecordSetBuilder();
        builder.addUnmappedFragment("read1");
        builder.addUnmappedFragment("read2");
        records = new ArrayList<>();
        builder.iterator().forEachRemaining(records::add);
    }

    @Test
    public void testPrefixRemoval() {
        records.forEach(record -> record.setAttribute("XT", "prefix_value"));
        underlyingIterator = records.iterator();
        BAMTagCleanupIterator iterator = new BAMTagCleanupIterator.Builder(underlyingIterator)
                .tag("XT")
                .prefixToRemove("prefix_")
                .build();

        SAMRecord result = iterator.next();
        assertEquals(result.getStringAttribute("XT"), "value");
    }

    @Test
    public void testSuffixRemoval() {
        records.forEach(record -> record.setAttribute("XT", "value_suffix"));
        underlyingIterator = records.iterator();
        BAMTagCleanupIterator iterator = new BAMTagCleanupIterator.Builder(underlyingIterator)
                .tag("XT")
                .suffixToRemove("_suffix")
                .build();

        SAMRecord result = iterator.next();
        assertEquals(result.getStringAttribute("XT"), "value");
    }

    @Test
    public void testPatternRemoval() {
        records.forEach(record -> record.setAttribute("XT", "value_to_remove"));
        underlyingIterator = records.iterator();
        BAMTagCleanupIterator iterator = new BAMTagCleanupIterator.Builder(underlyingIterator)
                .tag("XT")
                .patternToRemove(Pattern.compile("_to_remove"))
                .build();

        SAMRecord result = iterator.next();
        assertEquals(result.getStringAttribute("XT"), "value");
    }

    @Test
    public void testPrefixAndSuffixAddition() {
        records.forEach(record -> record.setAttribute("XT", "value"));
        underlyingIterator = records.iterator();
        BAMTagCleanupIterator iterator = new BAMTagCleanupIterator.Builder(underlyingIterator)
                .tag("XT")
                .prefixToAdd("prefix_")
                .suffixToAdd("_suffix")
                .build();

        SAMRecord result = iterator.next();
        assertEquals(result.getStringAttribute("XT"), "prefix_value_suffix");
    }

    @Test
    public void testCombinedTransformations() {
        records.forEach(record -> record.setAttribute("XT", "prefix_value_suffix"));
        underlyingIterator = records.iterator();
        BAMTagCleanupIterator iterator = new BAMTagCleanupIterator.Builder(underlyingIterator)
                .tag("XT")
                .prefixToRemove("prefix_")
                .suffixToRemove("_suffix")
                .prefixToAdd("new_prefix_")
                .suffixToAdd("_new_suffix")
                .build();

        SAMRecord result = iterator.next();
        assertEquals(result.getStringAttribute("XT"), "new_prefix_value_new_suffix");
    }
}