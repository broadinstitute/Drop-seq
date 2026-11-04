package org.broadinstitute.dropseqrna.utils;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class BatchingIterableTest {

    private static List<List<Integer>> batch(final int numItems, final int batchSize) {
        final List<Integer> items = IntStream.range(0, numItems).boxed().collect(Collectors.toList());
        final List<List<Integer>> result = new ArrayList<>();
        new BatchingIterable<>(items.iterator(), batchSize).forEach(result::add);
        return result;
    }

    @Test
    public void testLastBatchIsPartial() {
        Assert.assertEquals(batch(7, 3), Arrays.asList(
                Arrays.asList(0, 1, 2), Arrays.asList(3, 4, 5), Collections.singletonList(6)));
    }

    @Test
    public void testExactMultipleHasNoEmptyBatch() {
        Assert.assertEquals(batch(6, 3), Arrays.asList(Arrays.asList(0, 1, 2), Arrays.asList(3, 4, 5)));
    }

    @Test
    public void testFewerItemsThanBatchSize() {
        Assert.assertEquals(batch(2, 10), Collections.singletonList(Arrays.asList(0, 1)));
    }

    @Test
    public void testBatchSizeOne() {
        Assert.assertEquals(batch(3, 1), Arrays.asList(
                Collections.singletonList(0), Collections.singletonList(1), Collections.singletonList(2)));
    }

    @Test
    public void testEmptySourceGivesNoBatches() {
        Assert.assertTrue(batch(0, 3).isEmpty());
        final Iterator<List<Integer>> it = new BatchingIterable<>(Collections.<Integer>emptyIterator(), 3).iterator();
        Assert.assertFalse(it.hasNext());
        Assert.assertThrows(NoSuchElementException.class, it::next);
    }

    @Test
    public void testOrderIsPreserved() {
        final List<Integer> flattened = batch(1003, 100).stream().flatMap(List::stream).collect(Collectors.toList());
        Assert.assertEquals(flattened, IntStream.range(0, 1003).boxed().collect(Collectors.toList()));
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testBatchSizeMustBePositive() {
        new BatchingIterable<>(Collections.<Integer>emptyIterator(), 0);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testOnlyIterableOnce() {
        final BatchingIterable<Integer> batches = new BatchingIterable<>(Collections.<Integer>emptyIterator(), 3);
        batches.iterator();
        batches.iterator();
    }
}
