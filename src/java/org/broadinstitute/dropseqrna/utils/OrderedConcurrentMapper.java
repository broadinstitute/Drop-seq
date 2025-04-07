/*
 * MIT License
 *
 * Copyright 2025 Broad Institute
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

import htsjdk.samtools.SAMFileWriterImpl;
import htsjdk.samtools.util.CloserUtil;

import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.PriorityBlockingQueue;
import java.util.concurrent.Semaphore;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.AtomicReference;
import java.util.function.BiFunction;
import java.util.function.Consumer;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.StreamSupport;

/**
 * Concurrently iterate keeping the output in the order of the input.
 *
 * <p>One thread is used for consuming items from the queue and adding them to the final result.
 * The remaining threads are used for producing items from the iterable and adding them to the consumer's queue.
 * Does not allow more that a maximum number of items in the consumer queue at any time.
 * When done iterating, closes the iterable and any per thread instances using CloserUtil.
 *
 * <p>Pseudo code:
 * <pre>
 * StreamSupport
 *     .stream(iterable.spliterator(), true)
 *     .onClose(() -> CloserUtil.close(perThreadInstances))
 *     .onClose(() -> CloserUtil.close(iterable))
 *     // map with threads = (numThreads - 1)
 *     .map(input -> mapper(input, perThreadInstance))
 *     // consume with 1 dedicated thread, max items in the consumer's queue = maxItems
 *     .forEachOrdered(consumer);
 * </pre>
 */
public interface OrderedConcurrentMapper {

    /**
     * Run the iteration and wait for it to finish successfully or with a RuntimeException.
     */
    void execute();

    /**
     * Builder class for OrderedConcurrentMapper.
     *
     * @param <Input>     The type of the input being iterated and mapped
     * @param <PerThread> The type of the per-thread instance used during mapping
     * @param <Output>    The type returned by the mapper and then output
     */
    class Builder<Input, PerThread, Output> {
        private Iterable<Input> iterable;
        private int numThreads = 1;
        private int maxItems = SAMFileWriterImpl.getDefaultMaxRecordsInRam();
        private Supplier<PerThread> perThreadFactory;
        private BiFunction<Input, PerThread, Output> mapper;
        private Consumer<Output> consumer;

        /**
         * Sets the source iterable to be processed by the mapper.
         *
         * @param iterable The iterable containing input elements to process
         * @return This builder instance for method chaining
         */
        public Builder<Input, PerThread, Output> withIterable(final Iterable<Input> iterable) {
            this.iterable = iterable;
            return this;
        }

        /**
         * Sets the number of threads to use for parallel processing.
         * Must be at least 2 for parallel execution, as one thread is reserved for
         * consuming results while the others produce them.
         *
         * <p>Set to <1 to use the number of available processors.
         *
         * @param numThreads The number of threads to use (minimum 2 for parallel execution)
         * @return This builder instance for method chaining
         */
        public Builder<Input, PerThread, Output> withNumThreads(final int numThreads) {
            this.numThreads = numThreads;
            return this;
        }

        /**
         * Sets the maximum number of items that can be in the processing queue at once.
         * Higher values may use more memory but can improve throughput by reducing thread contention.
         *
         * @param maxItems The maximum number of items in the processing queue
         * @return This builder instance for method chaining
         */
        public Builder<Input, PerThread, Output> withMaxItems(final int maxItems) {
            this.maxItems = maxItems;
            return this;
        }

        /**
         * Sets the factory for creating per-thread context objects.
         * Each worker thread will receive its own instance of the context object.
         *
         * @param factory A supplier that creates new per-thread context instances
         * @return This builder instance for method chaining
         */
        public Builder<Input, PerThread, Output> withPerThreadFactory(final Supplier<PerThread> factory) {
            this.perThreadFactory = factory;
            return this;
        }

        /**
         * Sets the mapping function that transforms input items into output items.
         * The mapping function receives both the input item and the thread-local context object.
         *
         * @param mapper A function that maps from input to output using the per-thread context
         * @return This builder instance for method chaining
         */
        public Builder<Input, PerThread, Output> withMapper(final BiFunction<Input, PerThread, Output> mapper) {
            this.mapper = mapper;
            return this;
        }

        /**
         * Sets the mapping function that transforms input items into output items.
         * The mapping function receives both the input item and the thread-local context object.
         *
         * @param mapper A function that maps from input to without using the per-thread context
         * @return This builder instance for method chaining
         */
        public Builder<Input, PerThread, Output> withMapper(final Function<Input, Output> mapper) {
            if (mapper == null) {
                this.mapper = null;
            } else {
                this.mapper = (input, unused) -> mapper.apply(input);
            }
            return this;
        }

        /**
         * Sets the consumer that will receive the output items in the original iteration order.
         *
         * @param consumer A consumer that accepts output items
         * @return This builder instance for method chaining
         */
        public Builder<Input, PerThread, Output> withConsumer(final Consumer<Output> consumer) {
            this.consumer = consumer;
            return this;
        }

        /**
         * Validates the configuration, creates a MapperConfig, and builds an OrderedConcurrentMapper
         * instance that can be executed later by calling loopAndWait().
         *
         * @return A configured OrderedConcurrentMapper instance
         * @throws IllegalStateException if any required configuration is missing
         */
        public OrderedConcurrentMapper build() {
            if (iterable == null) throw new IllegalStateException("Iterable must be provided");
            if (mapper == null) throw new IllegalStateException("Mapper function must be provided");
            if (consumer == null) throw new IllegalStateException("Consumer must be provided");

            if (perThreadFactory == null) {
                perThreadFactory = () -> null; // Default to a no-op factory if none provided
            }

            if (numThreads < 1) {
                numThreads = Runtime.getRuntime().availableProcessors();
            }

            if (numThreads == 1) {
                return new SingleThreaded<>(iterable, perThreadFactory, mapper, consumer);
            } else {
                return new MultiThreaded<>(iterable, numThreads, maxItems, perThreadFactory, mapper, consumer);
            }
        }
    }

    /**
     * Single threaded implementation of the OrderedConcurrentMapper.
     *
     * @param <Input>     The type of the input being iterated and mapped
     * @param <PerThread> The type of the per-thread instance used during mapping
     * @param <Output>    The type returned by the mapper and then output
     */
    class SingleThreaded<Input, PerThread, Output> implements OrderedConcurrentMapper {
        /**
         * The iterable to process.
         */
        private final Iterable<Input> iterable;

        /**
         * Create a new per-thread instance. One is passed in to the mapper.
         */
        private final Supplier<PerThread> perThreadFactory;

        /**
         * Maps from Input to Output. The per-thread instance is passed in.
         */
        private final BiFunction<Input, PerThread, Output> mapper;

        /**
         * Adds the output to the final result. This is called after all items have been mapped.
         */
        private final Consumer<Output> consumer;

        public SingleThreaded(
                final Iterable<Input> iterable,
                final Supplier<PerThread> perThreadFactory,
                final BiFunction<Input, PerThread, Output> mapper,
                final Consumer<Output> consumer
        ) {
            this.iterable = iterable;
            this.perThreadFactory = perThreadFactory;
            this.mapper = mapper;
            this.consumer = consumer;
        }

        @Override
        public void execute() {
            final PerThread perThreadInstance = this.perThreadFactory.get();
            StreamSupport
                    .stream(this.iterable.spliterator(), false)
                    .onClose(() -> CloserUtil.close(perThreadInstance))
                    .onClose(() -> CloserUtil.close(this.iterable))
                    .map(input -> this.mapper.apply(input, perThreadInstance))
                    .forEachOrdered(this.consumer);
        }
    }

    /**
     * Multi-threaded implementation of the OrderedConcurrentMapper.
     *
     * @param <Input>     The type of the input being iterated and mapped
     * @param <PerThread> The type of the per-thread instance used during mapping
     * @param <Output>    The type returned by the mapper and then output
     */
    class MultiThreaded<Input, PerThread, Output> implements OrderedConcurrentMapper {

        /**
         * The iterable to process in parallel.
         */
        private final Iterable<Input> iterable;

        /**
         * Maps from Input to Output. One per-thread instance is passed in.
         */
        private final BiFunction<Input, PerThread, Output> mapper;

        /**
         * Adds the output to the final result. This is called after all items have been mapped.
         */
        private final Consumer<Output> consumer;

        /**
         * Number of threads to use for processing.
         */
        private final int numThreads;

        /**
         * ThreadLocal to keep track of the per-thread instance.
         */
        private final ThreadLocal<PerThread> perThreadLocal;

        /**
         * Concurrent collection to keep track of the per-thread instances.
         * Since the instances may not be suitable for storing in a Set, using ConcurrentLinkedQueue.
         */
        private final ConcurrentLinkedQueue<PerThread> perThreadInstances = new ConcurrentLinkedQueue<>();

        /**
         * Counter to keep track of the order of items produced from the iterator.
         */
        private final AtomicLong producerCounter = new AtomicLong();

        /**
         * Semaphore to control the number of items in the queue.
         */
        private final Semaphore consumerSemaphore;

        /**
         * Priority queue to keep track of the order of the outputs.
         */
        private final PriorityBlockingQueue<OrderedItem<Output>> consumerQueue;

        /**
         * Keep track of any exception that occurred.
         */
        private final AtomicReference<RuntimeException> exceptionRef = new AtomicReference<>();

        /**
         * Keep track if the loop has already been called.
         */
        private final AtomicBoolean alreadyCalled = new AtomicBoolean();

        public MultiThreaded(
                final Iterable<Input> iterable,
                final int numThreads,
                final int maxItems,
                final Supplier<PerThread> perThreadFactory,
                final BiFunction<Input, PerThread, Output> mapper,
                final Consumer<Output> consumer
        ) {
            if (numThreads < 2) {
                throw new IllegalArgumentException("numThreads must be at least 2: " + numThreads);
            }
            this.iterable = iterable;
            this.mapper = mapper;
            this.consumer = consumer;
            this.consumerSemaphore = new Semaphore(maxItems);
            this.consumerQueue = new PriorityBlockingQueue<>(maxItems);
            this.numThreads = numThreads;
            this.perThreadLocal = ThreadLocal.withInitial(() -> {
                final PerThread instance = tryRun(perThreadFactory::get);
                if (instance != null) {
                    this.perThreadInstances.add(instance);
                }
                return instance;
            });
        }

        @Override
        public void execute() {
            if (this.alreadyCalled.getAndSet(true)) {
                throw new IllegalStateException("loopAndWait() can only be called once.");
            }
            try (final ExecutorService producingExecutor =
                         Executors.newFixedThreadPool(this.numThreads - 1, new NamedThreadFactory("Producer"));
                 final ExecutorService consumerExecutor =
                         Executors.newSingleThreadExecutor(new NamedThreadFactory("Consumer"))) {

                consume(consumerExecutor, producingExecutor);
                produce(producingExecutor);

                producingExecutor.shutdown();
                consumerExecutor.shutdown();
            } finally {
                CloserUtil.close(this.perThreadInstances.stream().toList());
                CloserUtil.close(this.iterable);
            }
        }

        /**
         * Consume items from the queue and add them to the final result.
         * The producingExecutor is passed in so it may be checked for its termination.
         *
         * @param consumerExecutor  The executor to use for consuming items.
         * @param producingExecutor The executor to use for producing items.
         */
        private void consume(final ExecutorService consumerExecutor,
                             final ExecutorService producingExecutor) {
            consumerExecutor.submit(() -> {
                long expectedIndex = 0;
                while (true) {
                    checkException();

                    // Get the next item from the queue
                    final OrderedItem<Output> entry = tryRun(() -> this.consumerQueue.poll(1, TimeUnit.SECONDS));

                    if (entry == null) {
                        if (producingExecutor.isTerminated()) {
                            // All tasks are done and the queue is empty
                            break;
                        } else {
                            continue;
                        }
                    }

                    if (entry.index() == expectedIndex) {
                        tryRun(() -> this.consumer.accept(entry.item()));
                        expectedIndex++;
                        this.consumerSemaphore.release();
                    } else {
                        // Reinsert the item back into the queue
                        // There should still be space in the queue since we didn't release the semaphore yet
                        tryRun(() -> this.consumerQueue.put(entry));
                    }
                }
            });
        }

        /**
         * Produce items from the iterable and add them to the queue.
         * Closes the iterator when done.
         *
         * @param producingExecutor The executor to use for producing items.
         */
        private void produce(final ExecutorService producingExecutor) {
            for (final Input input : this.iterable) {
                final long itemIndex = this.producerCounter.getAndIncrement();

                // Block until space is available in the queue
                while (true) {
                    checkException();
                    final boolean available = tryRun(() -> this.consumerSemaphore.tryAcquire(1, TimeUnit.SECONDS));
                    if (available) {
                        break;
                    }
                }

                producingExecutor.submit(() -> {
                    checkException();

                    // Map the input to output
                    final Output output = tryRun(() -> this.mapper.apply(input, this.perThreadLocal.get()));

                    // Add the output to the queue
                    tryRun(() -> this.consumerQueue.put(new OrderedItem<>(itemIndex, output)));
                });
            }
        }

        /**
         * Run a runnable that can throw InterruptedException.
         *
         * @param runnable The runnable to run.
         */
        private void tryRun(final InterruptableRunnable runnable) {
            tryRun(
                    // Similar to wrapping with Executors.callable
                    // but our runnable can throw InterruptedException
                    () -> {
                        runnable.run();
                        return null;
                    }
            );
        }

        private <V> V tryRun(final InterruptableCallable<V> callable) {
            try {
                return callable.call();
            } catch (final InterruptedException e) {
                Thread.currentThread().interrupt();
                final RuntimeException re = new RuntimeException(e);
                markException(re);
                throw re;
            } catch (final RuntimeException e) {
                markException(e);
                throw e;
            }
        }

        /**
         * Check if an exception occurred on any of the threads and throw it if it did.
         */
        private void checkException() {
            if (this.exceptionRef.get() != null) {
                throw this.exceptionRef.get();
            }
        }

        /**
         * Mark the exception that occurred.
         *
         * @param runtimeException The exception that occurred.
         */
        private void markException(final RuntimeException runtimeException) {
            this.exceptionRef.accumulateAndGet(
                    runtimeException,
                    (existingEx, newEx) -> existingEx != null ? existingEx : newEx
            );
        }

        /**
         * Interface to allow a Runnable to throw InterruptedException.
         */
        private interface InterruptableRunnable {
            void run() throws InterruptedException;
        }

        /**
         * Interface to allow a Callable to throw InterruptedException.
         */
        private interface InterruptableCallable<V> {
            V call() throws InterruptedException;
        }

        /**
         * Class to keep track of the order of the variants.
         *
         * @param index The order of the item.
         * @param item  The item to be ordered.
         */
        private record OrderedItem<Item>(long index, Item item) implements Comparable<OrderedItem<Item>> {
            @Override
            public int compareTo(final OrderedItem<Item> other) {
                return Long.compare(this.index, other.index);
            }
        }

        /**
         * Thread factory to create threads with a specific name prefix.
         */
        private static class NamedThreadFactory implements ThreadFactory {
            private final String namePrefix;
            private int threadNumber = 1;

            public NamedThreadFactory(final String namePrefix) {
                this.namePrefix = namePrefix;
            }

            @Override
            public Thread newThread(@SuppressWarnings("NullableProblems") final Runnable r) {
                final Thread t = new Thread(r);
                t.setName(namePrefix + "-thread-" + threadNumber++);
                return t;
            }
        }
    }
}
