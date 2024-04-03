/*
 * MIT License
 *
 * Copyright 2024 Broad Institute
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

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;

@SuppressWarnings("resource")
public class IteratorOfIteratorsTest {

    @Test
    public void testIteratorOfIterators() {
        final List<Integer> list1 = Arrays.asList(1, 2, 3);
        final List<Integer> list2 = Arrays.asList(4, 5, 6);
        final List<Integer> list3 = Arrays.asList(7, 8, 9);
        final List<Integer> list4 = Arrays.asList(10, 11, 12);
        final List<Integer> list5 = Arrays.asList(13, 14, 15);
        final List<Integer> list6 = Arrays.asList(16, 17, 18);
        final List<Integer> list7 = Arrays.asList(19, 20, 21);
        final List<Integer> list8 = Arrays.asList(22, 23, 24);
        final List<Integer> list9 = Arrays.asList(25, 26, 27);
        final List<Integer> list10 = Arrays.asList(28, 29, 30);

        final IteratorOfIterators<Integer> iteratorOfIterators = new IteratorOfIterators<>(
                Arrays.asList(
                        list1.iterator(),
                        list2.iterator(),
                        list3.iterator(),
                        list4.iterator(),
                        list5.iterator(),
                        list6.iterator(),
                        list7.iterator(),
                        list8.iterator(),
                        list9.iterator(),
                        list10.iterator()
                ).iterator());

        final List<Integer> expected = Arrays.asList(
                1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                21, 22, 23, 24, 25, 26, 27, 28, 29, 30
        );
        final List<Integer> actual = new ArrayList<>();
        while (iteratorOfIterators.hasNext()) {
            actual.add(iteratorOfIterators.next());
        }
        Assert.assertEquals(actual, expected);
    }

    @Test
    public void testIteratorOfIteratorsEmpty() {
        final IteratorOfIterators<Integer> iteratorOfIterators = new IteratorOfIterators<>(Collections.emptyIterator());
        Assert.assertFalse(iteratorOfIterators.hasNext());
    }

    @Test
    public void testIteratorOfIteratorsEmptyInner() {
        final IteratorOfIterators<Integer> iteratorOfIterators = new IteratorOfIterators<>(
                List.of(Collections.<Integer>emptyIterator()).iterator()
        );
        Assert.assertFalse(iteratorOfIterators.hasNext());
    }

    @Test
    public void testIteratorOfIteratorsEmptyInner2() {
        final IteratorOfIterators<Integer> iteratorOfIterators = new IteratorOfIterators<>(
                Arrays.asList(Collections.<Integer>emptyIterator(), Collections.<Integer>emptyIterator()).iterator()
        );
        Assert.assertFalse(iteratorOfIterators.hasNext());
    }

    @Test
    public void testIteratorOfClosableIterator() {
        final List<Integer> list1 = Arrays.asList(1, 2, 3);
        final List<Integer> list2 = Arrays.asList(4, 5, 6);
        final List<Integer> list3 = Arrays.asList(7, 8, 9);

        final ClosableTestIterator<Integer> closableIterator1 = new ClosableTestIterator<>(list1.iterator());
        final ClosableTestIterator<Integer> closableIterator2 = new ClosableTestIterator<>(list2.iterator());
        final ClosableTestIterator<Integer> closableIterator3 = new ClosableTestIterator<>(list3.iterator());

        final IteratorOfIterators<Integer> iteratorOfIterators = new IteratorOfIterators<>(
                Arrays.<Iterator<Integer>>asList(
                        closableIterator1,
                        closableIterator2,
                        closableIterator3
                ).iterator());

        final List<Integer> expected = Arrays.asList(
                1, 2, 3, 4, 5, 6, 7, 8, 9
        );
        final List<Integer> actual = new ArrayList<>();
        while (iteratorOfIterators.hasNext()) {
            actual.add(iteratorOfIterators.next());
        }
        Assert.assertEquals(actual, expected);
        Assert.assertTrue(closableIterator1.isClosed());
        Assert.assertTrue(closableIterator2.isClosed());
        Assert.assertTrue(closableIterator3.isClosed());
    }

    @Test
    public void testIteratorOfClosableIteratorEmpty() {
        final ClosableTestIterator<Integer> closableIterator1 = new ClosableTestIterator<>(Collections.emptyIterator());
        final ClosableTestIterator<Integer> closableIterator2 = new ClosableTestIterator<>(Collections.emptyIterator());
        final ClosableTestIterator<Integer> closableIterator3 = new ClosableTestIterator<>(Collections.emptyIterator());

        final IteratorOfIterators<Integer> iteratorOfIterators = new IteratorOfIterators<>(
                Arrays.<Iterator<Integer>>asList(
                        closableIterator1,
                        closableIterator2,
                        closableIterator3
                ).iterator());

        Assert.assertFalse(iteratorOfIterators.hasNext());
        Assert.assertTrue(closableIterator1.isClosed());
        Assert.assertTrue(closableIterator2.isClosed());
        Assert.assertTrue(closableIterator3.isClosed());
    }

    @Test
    public void testIteratorsNotClosedIfHasNext() {
        final List<Integer> list1 = Arrays.asList(1, 2, 3);
        final List<Integer> list2 = Arrays.asList(4, 5, 6);
        final List<Integer> list3 = Arrays.asList(7, 8, 9);

        final ClosableTestIterator<Integer> closableIterator1 = new ClosableTestIterator<>(list1.iterator());
        final ClosableTestIterator<Integer> closableIterator2 = new ClosableTestIterator<>(list2.iterator());
        final ClosableTestIterator<Integer> closableIterator3 = new ClosableTestIterator<>(list3.iterator());

        final IteratorOfIterators<Integer> iteratorOfIterators = new IteratorOfIterators<>(
                Arrays.<Iterator<Integer>>asList(
                        closableIterator1,
                        closableIterator2,
                        closableIterator3
                ).iterator());

        Assert.assertTrue(iteratorOfIterators.hasNext());
        Assert.assertFalse(closableIterator1.isClosed());
        Assert.assertFalse(closableIterator2.isClosed());
        Assert.assertFalse(closableIterator3.isClosed());
    }

    @Test
    public void testIteratorsPartialIterationNotClosed() {
        // Iterate over the first five elements of a nine element iterator of iterators
        final List<Integer> list1 = Arrays.asList(1, 2, 3);
        final List<Integer> list2 = Arrays.asList(4, 5, 6);
        final List<Integer> list3 = Arrays.asList(7, 8, 9);

        final ClosableTestIterator<Integer> closableIterator1 = new ClosableTestIterator<>(list1.iterator());
        final ClosableTestIterator<Integer> closableIterator2 = new ClosableTestIterator<>(list2.iterator());
        final ClosableTestIterator<Integer> closableIterator3 = new ClosableTestIterator<>(list3.iterator());

        final IteratorOfIterators<Integer> iteratorOfIterators = new IteratorOfIterators<>(
                Arrays.<Iterator<Integer>>asList(
                        closableIterator1,
                        closableIterator2,
                        closableIterator3
                ).iterator());

        for (int i = 0; i < 5; i++) {
            iteratorOfIterators.next();
        }

        Assert.assertTrue(iteratorOfIterators.hasNext());
        Assert.assertEquals(iteratorOfIterators.next(), 6);
        Assert.assertTrue(closableIterator1.isClosed());
        Assert.assertFalse(closableIterator2.isClosed());
        Assert.assertFalse(closableIterator3.isClosed());
    }

    @Test
    public void testIteratorsPartialIterationExplicitClose() {
        // Iterate over the first five elements of a nine element iterator of iterators
        final List<Integer> list1 = Arrays.asList(1, 2, 3);
        final List<Integer> list2 = Arrays.asList(4, 5, 6);
        final List<Integer> list3 = Arrays.asList(7, 8, 9);

        final ClosableTestIterator<Integer> closableIterator1 = new ClosableTestIterator<>(list1.iterator());
        final ClosableTestIterator<Integer> closableIterator2 = new ClosableTestIterator<>(list2.iterator());
        final ClosableTestIterator<Integer> closableIterator3 = new ClosableTestIterator<>(list3.iterator());

        final IteratorOfIterators<Integer> iteratorOfIterators = new IteratorOfIterators<>(
                Arrays.<Iterator<Integer>>asList(
                        closableIterator1,
                        closableIterator2,
                        closableIterator3
                ).iterator());

        for (int i = 0; i < 5; i++) {
            iteratorOfIterators.next();
        }

        iteratorOfIterators.close();
        Assert.assertTrue(closableIterator1.isClosed());
        Assert.assertTrue(closableIterator2.isClosed());
        Assert.assertTrue(closableIterator3.isClosed());
    }

    @Test
    public void testIteratorsPartialIterationExplicitCloseNoop() {
        // Iterate over the first five elements of a nine element iterator of iterators
        final List<Integer> list1 = Arrays.asList(1, 2, 3);
        final List<Integer> list2 = Arrays.asList(4, 5, 6);
        final List<Integer> list3 = Arrays.asList(7, 8, 9);

        final IteratorOfIterators<Integer> iteratorOfIterators = new IteratorOfIterators<>(
                Arrays.asList(
                        list1.iterator(),
                        list2.iterator(),
                        list3.iterator()
                ).iterator());

        for (int i = 0; i < 5; i++) {
            iteratorOfIterators.next();
        }

        // Closing the outer iterator should not matter if the inner iterators are not closable
        iteratorOfIterators.close();
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testNextWithoutHasNext() {
        final List<Integer> list1 = Arrays.asList(1, 2, 3);
        final List<Integer> list2 = Arrays.asList(4, 5, 6);
        final List<Integer> list3 = Arrays.asList(7, 8, 9);

        final ClosableTestIterator<Integer> closableIterator1 = new ClosableTestIterator<>(list1.iterator());
        final ClosableTestIterator<Integer> closableIterator2 = new ClosableTestIterator<>(list2.iterator());
        final ClosableTestIterator<Integer> closableIterator3 = new ClosableTestIterator<>(list3.iterator());

        final IteratorOfIterators<Integer> iteratorOfIterators = new IteratorOfIterators<>(
                Arrays.<Iterator<Integer>>asList(
                        closableIterator1,
                        closableIterator2,
                        closableIterator3
                ).iterator());

        // Iterate over all elements
        while (iteratorOfIterators.hasNext()) {
            iteratorOfIterators.next();
        }

        // This should throw an exception
        iteratorOfIterators.next();
    }

    @Test
    public void testFromIterator() {
        final List<Integer> list1 = Arrays.asList(1, 2, 3);
        final List<Integer> list2 = Arrays.asList(4, 5, 6);
        final List<Integer> list3 = Arrays.asList(7, 8, 9);

        // Test the method fromIterator
        final IteratorOfIterators<Integer> iteratorOfIterators = IteratorOfIterators.fromIterator(
                Arrays.asList(list1, list2, list3).iterator(),
                List::iterator
        );

        final List<Integer> expected = Arrays.asList(1, 2, 3, 4, 5, 6, 7, 8, 9);
        final List<Integer> actual = new ArrayList<>();
        while (iteratorOfIterators.hasNext()) {
            actual.add(iteratorOfIterators.next());
        }
        Assert.assertEquals(actual, expected);
    }

    @Test
    public void testFromIteratorEmpty() {
        // Test the method fromIterator with an empty outer iterator
        final IteratorOfIterators<Integer> iteratorOfIterators = IteratorOfIterators.fromIterator(
                Collections.<List<Integer>>emptyIterator(),
                List::iterator
        );

        Assert.assertFalse(iteratorOfIterators.hasNext());
    }

    @Test
    public void testFromIteratorLazy() {

        final AtomicInteger counter = new AtomicInteger(0);

        // A class that allows us to check if the inner iterators are created lazily
        class LazyTestIterator<T> implements Iterator<T> {
            private final Iterator<T> internal;

            public LazyTestIterator(final Iterator<T> internal) {
                counter.incrementAndGet();
                this.internal = internal;
            }

            @Override
            public boolean hasNext() {
                return internal.hasNext();
            }

            @Override
            public T next() {
                return internal.next();
            }
        }

        final List<Integer> list1 = Arrays.asList(1, 2, 3);
        final List<Integer> list2 = Arrays.asList(4, 5, 6);
        final List<Integer> list3 = Arrays.asList(7, 8, 9);

        final IteratorOfIterators<Integer> iteratorOfIterators = IteratorOfIterators.fromIterator(
                Arrays.asList(list1.iterator(), list2.iterator(), list3.iterator()).iterator(),
                LazyTestIterator::new
        );

        // We should be able to iterate all the way to the sixth element without creating the third iterator.
        for (int i = 0; i < 6; i++) {
            iteratorOfIterators.next();
        }

        Assert.assertEquals(counter.get(), 2);
    }
}
