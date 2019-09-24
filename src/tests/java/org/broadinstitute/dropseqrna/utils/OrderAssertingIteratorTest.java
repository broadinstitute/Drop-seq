package org.broadinstitute.dropseqrna.utils;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

public class OrderAssertingIteratorTest {

    @Test
    public void testNext() {
        List<Integer> xs = Arrays.asList(1,2,3,4,5,6,7);
        Comparator<Integer> comparator = new Comparator<Integer>() {
            @Override
            public int compare(Integer o1, Integer o2) {
                return o1 - o2;
            }
        };
        OrderAssertingIterator<Integer> xiter = new OrderAssertingIterator<>(xs.iterator(), comparator);
        for (int x : xiter) { Assert.assertTrue(true); }
        xiter.close();

        xs = Arrays.asList(3,2,1);
        xiter = new OrderAssertingIterator<>(xs.iterator(), comparator);
        try {
            for (int x : xiter) { Assert.assertTrue(true); }
            Assert.fail("Did not throw proper assertion");
        } catch (IllegalStateException e) {
            Assert.assertTrue(true);
        }
        xiter.close();

        xs = Arrays.asList(1,2,3,4,7,5,6);
        xiter = new OrderAssertingIterator<>(xs.iterator(), comparator);
        try {
            for (int x : xiter) { Assert.assertTrue(true); }
            Assert.fail("Did not throw proper assertion");
        } catch (IllegalStateException e) {
            Assert.assertTrue(true);
        }
        xiter.close();

        xs = Arrays.asList(7,1,2,3,4,5,6);
        xiter = new OrderAssertingIterator<>(xs.iterator(), comparator);
        try {
            for (int x : xiter) { Assert.assertTrue(true); }
            Assert.fail("Did not throw proper assertion");
        } catch (IllegalStateException e) {
            Assert.assertTrue(true);
        }
        xiter.close();
    }
}