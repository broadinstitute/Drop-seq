/*
 * MIT License
 *
 * Copyright 2026 Broad Institute
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
package org.broadinstitute.dropseqrna.utils.batchwriter;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicInteger;

public class BatchSinkGroupTest {

    /** What one sink received, and whether its batches ever overlapped. */
    private static class Recorder {
        final List<Integer> items = Collections.synchronizedList(new ArrayList<>());
        final List<String> events = Collections.synchronizedList(new ArrayList<>());
        final AtomicInteger active = new AtomicInteger();
        final AtomicBoolean overlapped = new AtomicBoolean();

        void write(final List<Integer> batch) {
            if (active.incrementAndGet() > 1) {
                overlapped.set(true);
            }
            items.addAll(batch);
            events.add("batch");
            active.decrementAndGet();
        }
    }

    @Test(timeOut = 60000)
    public void testOrderIsKeptPerSinkWithFewerThreadsThanSinks() {
        final int numSinks = 8;
        final int numItems = 200000;
        final BatchSinkGroup group = new BatchSinkGroup(3, 7, 2);
        final Recorder[] recorders = new Recorder[numSinks];
        final List<BatchSink<Integer>> sinks = new ArrayList<>();
        for (int i = 0; i < numSinks; i++) {
            final Recorder recorder = new Recorder();
            recorders[i] = recorder;
            sinks.add(group.open(recorder::write, () -> recorder.events.add("closed")));
        }
        for (int i = 0; i < numItems; i++) {
            sinks.get(i % numSinks).add(i);
        }
        group.close();
        for (int s = 0; s < numSinks; s++) {
            final List<Integer> expected = new ArrayList<>();
            for (int i = s; i < numItems; i += numSinks) {
                expected.add(i);
            }
            Assert.assertEquals(recorders[s].items, expected, "sink " + s);
            Assert.assertFalse(recorders[s].overlapped.get(), "batches of sink " + s + " overlapped");
        }
    }

    @Test(timeOut = 60000)
    public void testPartialBatchIsFlushedOnCloseBeforeOnClose() {
        final BatchSinkGroup group = new BatchSinkGroup(2, 100, 2);
        final Recorder recorder = new Recorder();
        final BatchSink<Integer> sink = group.open(recorder::write, () -> recorder.events.add("closed"));
        for (int i = 0; i < 5; i++) {
            sink.add(i);
        }
        Assert.assertTrue(recorder.items.isEmpty(), "a batch that is not full is not handed over");
        sink.close();
        Assert.assertEquals(recorder.items, List.of(0, 1, 2, 3, 4));
        Assert.assertEquals(recorder.events, List.of("batch", "closed"));
        group.close();
    }

    @Test(timeOut = 60000)
    public void testFailureInBatchWriterSurfacesAndDoesNotHang() {
        final BatchSinkGroup group = new BatchSinkGroup(2, 1, 2);
        final AtomicInteger batches = new AtomicInteger();
        final AtomicBoolean closed = new AtomicBoolean();
        final BatchSink<Integer> sink = group.open(b -> {
            if (batches.incrementAndGet() == 3) {
                throw new IllegalArgumentException("boom");
            }
        }, () -> closed.set(true));
        final RuntimeException e = Assert.expectThrows(RuntimeException.class, () -> {
            for (int i = 0; i < 100000; i++) {
                sink.add(i);
            }
            sink.close();
        });
        Assert.assertTrue(e.getMessage().contains("boom"), e.getMessage());
        Assert.assertEquals(batches.get(), 3, "nothing is processed after the failure");
        Assert.assertFalse(closed.get(), "onClose is not run after a failure");
    }

    @Test(timeOut = 60000)
    public void testProducerWaitsWhenBatchesAreNotBeingProcessed() throws InterruptedException {
        final int maxInFlight = 2;
        final BatchSinkGroup group = new BatchSinkGroup(1, 1, maxInFlight);
        final CountDownLatch release = new CountDownLatch(1);
        final Recorder recorder = new Recorder();
        final BatchSink<Integer> sink = group.open(batch -> {
            try {
                release.await();
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }
            recorder.write(batch);
        }, () -> {});
        final AtomicInteger added = new AtomicInteger();
        final Thread producer = new Thread(() -> {
            for (int i = 0; i < 100; i++) {
                sink.add(i);
                added.incrementAndGet();
            }
        });
        producer.start();
        Thread.sleep(500);
        Assert.assertTrue(producer.isAlive(), "producer should be blocked");
        Assert.assertEquals(added.get(), maxInFlight, "only maxBatchesInFlight batches can be outstanding");
        release.countDown();
        producer.join();
        group.close();
        final List<Integer> expected = new ArrayList<>();
        for (int i = 0; i < 100; i++) {
            expected.add(i);
        }
        Assert.assertEquals(recorder.items, expected);
    }
}
