package org.broadinstitute.dropseqrna.utils;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.function.Function;
import java.util.stream.IntStream;

public class OrderedConcurrentMapperTest {
    @Test
    public void testExecuteSingleThread() {
        final List<Integer> input = IntStream.range(0, 10000).boxed().toList();
        final TestMapperBuilder mapperBuilder = new TestMapperBuilder(input, 1, 100);
        mapperBuilder.execute();
        Assert.assertEquals(mapperBuilder.getResults(), input);
    }

    @Test
    public void testExecuteAndWait() {
        final List<Integer> input = IntStream.range(0, 10000).boxed().toList();
        final List<Integer> expected = IntStream.range(1, 10001).boxed().toList();
        final TestMapperBuilder mapperBuilder = new TestMapperBuilder(input, 3, 100) {
            @Override
            public Integer map(final Integer integer, final Object unused) {
                return integer + 1;
            }
        };
        mapperBuilder.execute();
        Assert.assertEquals(mapperBuilder.getResults(), expected);
    }

    @Test
    public void testExecuteMultiThread() {
        final List<Integer> input = IntStream.range(0, 10000).boxed().toList();
        final TestMapperBuilder mapperBuilder = new TestMapperBuilder(input, 3, 100);
        mapperBuilder.execute();
        Assert.assertEquals(mapperBuilder.getResults(), input);
    }

    @Test
    public void testExecuteNoPerThread() {
        final List<Integer> input = IntStream.range(0, 10000).boxed().toList();
        final TestMapperBuilder mapperBuilder = new TestMapperBuilder(input, 3, 100);
        final OrderedConcurrentMapper.Builder<Integer, Object, Integer> builder = mapperBuilder.builder();
        builder.withPerThreadFactory(null)
                .withMapper(in -> mapperBuilder.map(in, null))
                .build()
                .execute();
        Assert.assertEquals(mapperBuilder.getResults(), input);
    }

    @Test
    public void testExecuteRandom() {
        final List<Integer> input = new Random(777).ints(10000).boxed().toList();
        final TestMapperBuilder mapperBuilder = new TestMapperBuilder(input, 3, 100);
        mapperBuilder.execute();
        Assert.assertEquals(mapperBuilder.getResults(), input);
    }

    @Test
    public void testSlowOutput() {
        final int maxItems = 3;
        final List<Integer> input = IntStream.range(0, 10000).boxed().toList();
        final AtomicInteger max = new AtomicInteger();
        final AtomicInteger cnt = new AtomicInteger();
        final TestMapperBuilder mapperBuilder = new TestMapperBuilder(input, maxItems * 2, maxItems) {
            @Override
            public Integer map(final Integer integer, final Object unused) {
                max.accumulateAndGet(cnt.incrementAndGet(), Math::max);
                try {
                    return super.map(integer, unused);
                } finally {
                    cnt.decrementAndGet();
                }
            }

            @Override
            public void add(final Integer integer) {
                try {
                    if (integer % 3000 == 0) {
                        Thread.sleep(100);
                    }
                } catch (final InterruptedException e) {
                    throw new RuntimeException(e);
                }
                super.add(integer);
            }
        };
        mapperBuilder.execute();
        Assert.assertEquals(mapperBuilder.getResults(), input);
        // Even with backup, we shouldn't have more than maxItems mapping at a time.
        Assert.assertTrue(max.get() > 0);
        Assert.assertTrue(max.get() <= maxItems);
        Assert.assertEquals(cnt.get(), 0);
    }

    @Test
    public void testPerThreadInstanceCount() {
        final int maxItems = 3;
        final int numProducerThreads = maxItems * 2;
        final List<Integer> input = IntStream.range(0, 10000).boxed().toList();
        final AtomicInteger cnt = new AtomicInteger();
        final TestMapperBuilder mapperBuilder = new TestMapperBuilder(input, numProducerThreads + 1, maxItems) {
            @Override
            public Object newPerThreadInstance() {
                cnt.incrementAndGet();
                return new Object();
            }
        };
        mapperBuilder.execute();
        Assert.assertEquals(mapperBuilder.getResults(), input);
        Assert.assertTrue(cnt.get() <= numProducerThreads);
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testExecuteTwice() {
        final List<Integer> input = IntStream.range(0, 10000).boxed().toList();
        final TestMapperBuilder mapperBuilder = new TestMapperBuilder(input, 2, 100);
        final OrderedConcurrentMapper mapper = mapperBuilder.build();
        mapper.execute();
        mapper.execute();
    }

    @Test(expectedExceptions = IllegalArgumentException.class)
    public void testExceptionOneMultiThread() {
        final List<Integer> input = IntStream.range(0, 10000).boxed().toList();
        final TestMapperBuilder mapperBuilder = new TestMapperBuilder(input, 1, 100);
        new OrderedConcurrentMapper.MultiThreaded<>(
                mapperBuilder.getIterable(),
                mapperBuilder.getNumThreads(),
                mapperBuilder.getMaxItems(),
                mapperBuilder::newPerThreadInstance,
                mapperBuilder::map,
                mapperBuilder::add
        );
    }

    @Test(expectedExceptions = ExpectedException.class)
    public void testNewPerThreadInstanceException() {
        final List<Integer> input = IntStream.range(0, 10000).boxed().toList();
        final TestMapperBuilder mapperBuilder = new TestMapperBuilder(input, 3, 100) {
            @Override
            public Object newPerThreadInstance() {
                throw new ExpectedException("testNewPerThreadInstanceException");
            }
        };
        mapperBuilder.execute();
    }

    @Test(expectedExceptions = ExpectedException.class)
    public void testMapException() {
        final List<Integer> input = IntStream.range(0, 10000).boxed().toList();
        final TestMapperBuilder mapperBuilder = new TestMapperBuilder(input, 3, 100) {
            @Override
            public Integer map(final Integer integer, final Object unused) {
                if (integer == 3000) {
                    throw new ExpectedException("testMapException");
                }
                return integer;
            }
        };
        mapperBuilder.execute();
    }

    @Test(expectedExceptions = ExpectedException.class)
    public void testAddException() {
        final List<Integer> input = IntStream.range(0, 10000).boxed().toList();
        final TestMapperBuilder mapperBuilder = new TestMapperBuilder(input, 3, 100) {
            @Override
            public void add(final Integer integer) {
                if (integer == 3000) {
                    throw new ExpectedException("testAddException");
                }
                super.add(integer);
            }
        };
        mapperBuilder.execute();
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testBuilderNoIterator() {
        final List<Integer> input = IntStream.range(0, 10000).boxed().toList();
        final TestMapperBuilder mapperBuilder = new TestMapperBuilder(input, 3, 100);
        final OrderedConcurrentMapper.Builder<Integer, Object, Integer> builder = mapperBuilder.builder();
        builder.withIterable(null);
        builder.build();
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testBuilderNoMapper() {
        final List<Integer> input = IntStream.range(0, 10000).boxed().toList();
        final TestMapperBuilder mapperBuilder = new TestMapperBuilder(input, 3, 100);
        final OrderedConcurrentMapper.Builder<Integer, Object, Integer> builder = mapperBuilder.builder();
        builder.withMapper((Function<Integer, Integer>) null);
        builder.build();
    }

    @Test(expectedExceptions = IllegalStateException.class)
    public void testBuilderNoConsumer() {
        final List<Integer> input = IntStream.range(0, 10000).boxed().toList();
        final TestMapperBuilder mapperBuilder = new TestMapperBuilder(input, 3, 100);
        final OrderedConcurrentMapper.Builder<Integer, Object, Integer> builder = mapperBuilder.builder();
        builder.withConsumer(null);
        builder.build();
    }

    private static class TestMapperBuilder {
        private final List<Integer> results = new ArrayList<>();
        private final List<Integer> input;
        private final int numThreads;
        private final int maxItems;

        public TestMapperBuilder(final List<Integer> input, final int numThreads, final int maxItems) {
            this.input = input;
            this.numThreads = numThreads;
            this.maxItems = maxItems;
        }

        public List<Integer> getResults() {
            return results;
        }

        public Iterable<Integer> getIterable() {
            return input;
        }

        public int getNumThreads() {
            return numThreads;
        }

        public int getMaxItems() {
            return maxItems;
        }

        public Object newPerThreadInstance() {
            return null;
        }

        public Integer map(final Integer integer, final Object unused) {
            return integer;
        }

        public void add(final Integer integer) {
            results.add(integer);
        }

        public OrderedConcurrentMapper.Builder<Integer, Object, Integer> builder() {
            return new OrderedConcurrentMapper
                    .Builder<Integer, Object, Integer>()
                    .withNumThreads(this.getNumThreads())
                    .withMaxItems(this.getMaxItems())
                    .withIterable(this.getIterable())
                    .withPerThreadFactory(this::newPerThreadInstance)
                    .withMapper(this::map)
                    .withConsumer(this::add);
        }

        public OrderedConcurrentMapper build() {
            return builder().build();
        }

        public void execute() {
            builder().build().execute();
        }
    }

    private static class ExpectedException extends RuntimeException {
        public ExpectedException(final String message) {
            super("ExpectedException: " + message);
        }
    }
}
