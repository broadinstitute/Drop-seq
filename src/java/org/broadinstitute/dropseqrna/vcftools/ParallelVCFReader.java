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
package org.broadinstitute.dropseqrna.vcftools;

import htsjdk.samtools.util.BlockCompressedInputStream;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.tribble.readers.LineIteratorImpl;
import htsjdk.tribble.readers.PositionalBufferedStream;
import htsjdk.tribble.readers.SynchronousLineReader;
import htsjdk.variant.bcf2.BCF2Codec;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.LazyGenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.Closeable;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.io.UncheckedIOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;
import java.util.function.IntConsumer;

/**
 * Reads a VCF or BCF like {@link VCFFileReader}, but decodes the records on several threads.
 * The records come back in the same order, with the same content, as {@code VCFFileReader.iterator()}.
 *
 * <pre>{@code
 * try (ParallelVCFReader reader = new ParallelVCFReader(vcfPath, numThreads)) {
 *     VCFHeader header = reader.getFileHeader();
 *     for (VariantContext vc : reader) { ... }
 * }
 * }</pre>
 *
 * <p>One thread reads the file and splits it into batches of records.  Worker threads decode the batches, and
 * the iterator returns them in file order.  With one thread, this just uses {@link VCFFileReader}.
 *
 * <p><b>Differences from {@link VCFFileReader}:</b> this class only reads the file from start to end; it has no
 * {@code query()} to jump to a region.  Iterating never needs the index, with either reader.  Like {@code VCFFileReader},
 * the constructors fail if the index is missing unless you pass {@code requireIndex = false}, but with more than one
 * thread that is only an existence check.  With one thread this simply wraps {@link VCFFileReader}.
 *
 * <p>Supported files: VCF (plain or gzipped) and uncompressed BCF (named {@code *.bcf}).
 * BGZF-compressed BCF is not supported until htsjdk 6.
 *
 * <p>With more than one thread, each record's genotypes are fully decoded before it is returned.  htsjdk normally
 * decodes genotypes lazily, but that isn't safe across threads.
 *
 * <p>More threads help only while decoding is the slowest step.  Text VCFs, which are slow to decode, benefit
 * the most.  BCF decodes much faster, so the caller's own work on each record soon becomes the limit.
 *
 * <p>Tuning (the defaults are usually fine):
 * <ul>
 *   <li>{@code batchSize}: records per batch (default {@value #DEFAULT_BATCH_SIZE}).  Very large batches are slower.</li>
 *   <li>{@code queueDepth}: batches decoded ahead of the caller (default 2 × numThreads).  Use at least numThreads.</li>
 *   <li>Memory: up to about {@code (queueDepth + 2) × batchSize} decoded records are held at once.</li>
 * </ul>
 *
 * <p>Errors from reading or decoding are thrown by the iterator.  Always close the reader, even if you stop early.
 */
public class ParallelVCFReader implements Closeable, Iterable<VariantContext> {

	private static final Log log = Log.getInstance(ParallelVCFReader.class);
	private static final Future<List<VariantContext>> END = CompletableFuture.completedFuture(Collections.emptyList());
	private static final AtomicInteger INSTANCE_COUNTER = new AtomicInteger();
	private static final int INPUT_BUFFER_SIZE = 1 << 20;
	// How long close() waits for each of the reader thread and the workers to stop.  They are daemon threads, so one that
	// doesn't stop in time (a stuck read, say) is reported and abandoned rather than blocking close() or JVM exit.
	private static final long CLOSE_WAIT_SECONDS = 10;

	/** Default number of records per batch. */
	public static final int DEFAULT_BATCH_SIZE = 100;

	private final String threadNamePrefix;
	private final IntConsumer beforeBatchDecode;
	private final VCFHeader header;
	private final int batchSize;
	/** Non-null only when numThreads == 1, in which case every field below is null. */
	private final VCFFileReader serialReader;
	private final Format format;
	private final BlockingQueue<Future<List<VariantContext>>> queue;
	private final ExecutorService workers;
	private final Thread readerThread;
	private volatile boolean closed = false;
	private boolean iteratorCreated = false;

	// diagnostics
	private final AtomicLong consumerWaitNanos = new AtomicLong();
	private final AtomicLong readerBlockedNanos = new AtomicLong();

	/**
	 * Opens a VCF or BCF using the default tuning.  As with {@link VCFFileReader}, an index is required; use
	 * {@link #ParallelVCFReader(Path, boolean, int)} with {@code requireIndex = false} to read a file without one.
	 * @param vcf the VCF or BCF to read
	 * @param numThreads number of decoding threads; 1 means read with {@link VCFFileReader}
	 */
	public ParallelVCFReader(final Path vcf, final int numThreads) {
		this(vcf, true, numThreads);
	}

	/**
	 * Opens a VCF or BCF using the default tuning.
	 * @param vcf the VCF or BCF to read
	 * @param requireIndex if true, fail if the file has no index, as {@link VCFFileReader} does.
	 * Reading with more than one thread never uses the index; this is only an existence check.
	 * @param numThreads number of decoding threads; 1 means read with {@link VCFFileReader}
	 */
	public ParallelVCFReader(final Path vcf, final boolean requireIndex, final int numThreads) {
		this(vcf, requireIndex, numThreads, DEFAULT_BATCH_SIZE, 2 * numThreads);
	}

	/**
	 * Opens a VCF or BCF with explicit tuning.  The tuning is ignored when numThreads is 1.
	 * @param vcf the VCF or BCF to read
	 * @param requireIndex if true, fail if the file has no index, as {@link VCFFileReader} does.
	 * Reading with more than one thread never uses the index; this is only an existence check.
	 * @param numThreads number of decoding threads; 1 means read with {@link VCFFileReader}
	 * @param batchSize records per batch
	 * @param queueDepth batches decoded ahead of the caller; use at least numThreads
	 */
	public ParallelVCFReader(final Path vcf, final boolean requireIndex, final int numThreads, final int batchSize, final int queueDepth) {
		this(vcf, requireIndex, numThreads, batchSize, queueDepth, null);
	}

	/**
	 * Test hook: {@code beforeBatchDecode} is called on the worker thread with the batch's 0-based index
	 * just before the batch is decoded.  Lets tests force batches to finish out of order.
	 */
	ParallelVCFReader(final Path vcf, final boolean requireIndex, final int numThreads, final int batchSize, final int queueDepth,
					  final IntConsumer beforeBatchDecode) {
		if (numThreads < 1) throw new IllegalArgumentException("numThreads must be >= 1");
		if (batchSize < 1) throw new IllegalArgumentException("batchSize must be >= 1");
		if (queueDepth < 1) throw new IllegalArgumentException("queueDepth must be >= 1");
		this.threadNamePrefix = "ParallelVCFReader-" + INSTANCE_COUNTER.incrementAndGet() + "-";
		this.beforeBatchDecode = beforeBatchDecode;
		this.batchSize = batchSize;

		if (numThreads == 1) {
			// One thread: use htsjdk's reader directly.  No extra threads, and genotypes stay lazily decoded.
			this.serialReader = new VCFFileReader(vcf, requireIndex);
			this.header = serialReader.getFileHeader();
			this.format = null;
			this.queue = null;
			this.workers = null;
			this.readerThread = null;
			return;
		}
		this.serialReader = null;
		if (requireIndex) {
			// The parallel path never uses the index, but fail on a missing one the same way VCFFileReader does.
			new VCFFileReader(vcf, true).close();
		}
		try {
			this.format = VCFFileReader.isBCF(vcf) ? new BcfFormat(vcf) : new VcfFormat(vcf);
		} catch (IOException e) {
			throw new UncheckedIOException(e);
		}
		this.header = format.header();
		this.queue = new ArrayBlockingQueue<>(queueDepth);

		final AtomicInteger workerNum = new AtomicInteger();
		this.workers = Executors.newFixedThreadPool(numThreads, r -> {
			Thread t = new Thread(r, threadNamePrefix + "worker-" + workerNum.incrementAndGet());
			t.setDaemon(true);
			return t;
		});
		this.readerThread = new Thread(this::readLoop, threadNamePrefix + "reader");
		this.readerThread.setDaemon(true);
		log.info(String.format("Reading %s with %d decoding threads (batch size %d, queue depth %d)",
				vcf, numThreads, batchSize, queueDepth));
	}

	/** @return the file's header */
	public VCFHeader getFileHeader() {
		return header;
	}

	/**
	 * @return nanoseconds the caller has spent waiting for decoded records.  If this is near zero,
	 * more threads won't help.
	 */
	public long getConsumerWaitNanos() { return consumerWaitNanos.get(); }

	/**
	 * @return nanoseconds the reader thread has spent waiting for room in the queue.  If this is near zero,
	 * reading the file is the bottleneck.
	 */
	public long getReaderBlockedNanos() { return readerBlockedNanos.get(); }

	/** Prefix of the names of every thread this instance starts.  For tests. */
	String getThreadNamePrefix() { return threadNamePrefix; }

	private void readLoop() {
		try {
			int batchIndex = 0;
			Callable<List<VariantContext>> task;
			while (!closed && (task = format.nextBatch(batchSize)) != null)
				submit(task, batchIndex++);
			put(END);
		} catch (InterruptedException e) {
			// closed while blocked on a full queue
		} catch (Throwable t) {
			try { put(CompletableFuture.failedFuture(t)); } catch (InterruptedException ignored) { }
		}
	}

	private void submit(final Callable<List<VariantContext>> task, final int batchIndex) throws InterruptedException {
		put(workers.submit(() -> {
			if (beforeBatchDecode != null) beforeBatchDecode.accept(batchIndex);
			return task.call();
		}));
	}

	private void put(final Future<List<VariantContext>> f) throws InterruptedException {
		final long t0 = System.nanoTime();
		queue.put(f);
		readerBlockedNanos.addAndGet(System.nanoTime() - t0);
	}

	/** Forces the lazily-parsed genotypes to be decoded now, on the calling (worker) thread. */
	private static VariantContext fullyDecodeGenotypes(final VariantContext vc) {
		final GenotypesContext gc = vc.getGenotypes();
		if (gc instanceof LazyGenotypesContext lazy) lazy.decode();
		return vc;
	}

	/**
	 * Starts reading and returns the records in file order.  Can be called only once.
	 * @throws IllegalStateException if called a second time
	 */
	@Override
	public synchronized CloseableIterator<VariantContext> iterator() {
		if (iteratorCreated) throw new IllegalStateException("iterator() may only be called once");
		iteratorCreated = true;
		if (serialReader != null) return serialReader.iterator();
		readerThread.start();
		return new OrderedIterator();
	}

	/** Stops all threads and closes the file.  Safe to call before all records are read, or more than once. */
	@Override
	public void close() {
		if (closed) return;
		closed = true;
		if (serialReader != null) {
			serialReader.close();
			return;
		}
		readerThread.interrupt();
		workers.shutdownNow();
		Future<List<VariantContext>> f;
		while ((f = queue.poll()) != null) f.cancel(true);
		try {
			readerThread.join(TimeUnit.SECONDS.toMillis(CLOSE_WAIT_SECONDS));
			if (readerThread.isAlive())
				log.warn("Reader thread ", readerThread.getName(), " did not stop within ", CLOSE_WAIT_SECONDS, "s of close(); abandoning it.");
			if (!workers.awaitTermination(CLOSE_WAIT_SECONDS, TimeUnit.SECONDS))
				log.warn("Decoding threads did not stop within ", CLOSE_WAIT_SECONDS, "s of close(); abandoning them.");
		} catch (InterruptedException e) {
			Thread.currentThread().interrupt();
		}
		try {
			format.close();
		} catch (IOException e) {
			throw new UncheckedIOException(e);
		}
	}

	/**
	 * The format-specific parts: reading the header, and cutting the input into batches of raw records
	 * (on the reader thread) that a worker can decode independently with its own codec.
	 */
	private interface Format extends Closeable {
		VCFHeader header();

		/**
		 * Reads up to {@code batchSize} raw records and returns a task that decodes them, or null at end of input.
		 * Called only on the reader thread; the returned task runs on a worker thread.
		 */
		Callable<List<VariantContext>> nextBatch(int batchSize) throws IOException;
	}

	/** Text VCF, plain or gzipped/BGZF.  Raw records are lines. */
	private static final class VcfFormat implements Format {
		private final BufferedReader lineReader;
		private final List<String> headerLines = new ArrayList<>();
		private final VCFHeader header;
		private final ThreadLocal<VCFCodec> codecs = ThreadLocal.withInitial(this::newCodec);

		VcfFormat(final Path vcf) throws IOException {
			final InputStream in = IOUtil.isBlockCompressed(vcf)
					? new BlockCompressedInputStream(Files.newInputStream(vcf))
					: IOUtil.openFileForReading(vcf);
			// ISO-8859-1 maps each byte to one char, matching htsjdk's own line readers.
			this.lineReader = new BufferedReader(new InputStreamReader(in, StandardCharsets.ISO_8859_1), INPUT_BUFFER_SIZE);
			String line;
			while ((line = lineReader.readLine()) != null) {
				headerLines.add(line);
				if (line.startsWith("#CHROM")) break;
			}
			this.header = newCodec().getHeader();
		}

		/** Each worker parses its own copy of the header, so no header or codec state is shared across threads. */
		private VCFCodec newCodec() {
			final VCFCodec codec = new VCFCodec();
			codec.readActualHeader(new LineIteratorImpl(new SynchronousLineReader(new StringReader(String.join("\n", headerLines) + "\n"))));
			return codec;
		}

		@Override
		public VCFHeader header() { return header; }

		@Override
		public Callable<List<VariantContext>> nextBatch(final int batchSize) throws IOException {
			final List<String> lines = new ArrayList<>(batchSize);
			String line;
			while (lines.size() < batchSize && (line = lineReader.readLine()) != null)
				lines.add(line);
			if (lines.isEmpty()) return null;
			return () -> {
				final VCFCodec codec = codecs.get();
				final List<VariantContext> result = new ArrayList<>(lines.size());
				for (String l : lines) {
					final VariantContext vc = codec.decode(l);
					if (vc != null) result.add(fullyDecodeGenotypes(vc));
				}
				return result;
			};
		}

		@Override
		public void close() throws IOException { lineReader.close(); }
	}

	/**
	 * Uncompressed BCF (as Picard/htsjdk write it).  The header is the magic bytes, a little-endian uint32 length
	 * and the text header.  Each record is a little-endian uint32 l_shared, uint32 l_indiv, then l_shared + l_indiv
	 * bytes, so the reader thread can find record boundaries without decoding anything.
	 *
	 * <p>BGZF-compressed BCF (bcftools' default output) is rejected, as it is by htsjdk 5's {@link VCFFileReader}.
	 * htsjdk 6 reads it; supporting it here would only need a {@link BlockCompressedInputStream} under {@code in}.
	 */
	private static final class BcfFormat implements Format {
		private final InputStream in;
		private final byte[] headerBytes;
		private final VCFHeader header;
		private final ThreadLocal<BCF2Codec> codecs = ThreadLocal.withInitial(this::newCodec);
		private final byte[] lengths = new byte[8];

		BcfFormat(final Path bcf) throws IOException {
			if (IOUtil.isBlockCompressed(bcf))
				throw new IllegalArgumentException("BGZF-compressed BCF is not supported (it requires htsjdk 6): " + bcf);
			this.in = new BufferedInputStream(Files.newInputStream(bcf), INPUT_BUFFER_SIZE);
			final byte[] magicAndLength = in.readNBytes(9);
			if (magicAndLength.length < 9 || magicAndLength[0] != 'B' || magicAndLength[1] != 'C' || magicAndLength[2] != 'F')
				throw new IOException("Not a BCF file (bad magic bytes): " + bcf);
			final long textLength = Integer.toUnsignedLong(ByteBuffer.wrap(magicAndLength, 5, 4).order(ByteOrder.LITTLE_ENDIAN).getInt());
			if (textLength > Integer.MAX_VALUE - 9) throw new IOException("BCF header too large: " + bcf);
			final byte[] text = in.readNBytes((int) textLength);
			if (text.length < textLength) throw new EOFException("Truncated BCF header: " + bcf);
			this.headerBytes = new byte[9 + text.length];
			System.arraycopy(magicAndLength, 0, headerBytes, 0, 9);
			System.arraycopy(text, 0, headerBytes, 9, text.length);
			// Parsing the header here also checks that htsjdk supports this BCF version.
			this.header = (VCFHeader) newCodecAndHeader().header;
		}

		private record CodecAndHeader(BCF2Codec codec, Object header) { }

		private CodecAndHeader newCodecAndHeader() {
			final BCF2Codec codec = new BCF2Codec();
			final Object h = codec.readHeader(new PositionalBufferedStream(new ByteArrayInputStream(headerBytes))).getHeaderValue();
			return new CodecAndHeader(codec, h);
		}

		/** Each worker parses its own copy of the header, so no header or codec state is shared across threads. */
		private BCF2Codec newCodec() {
			return newCodecAndHeader().codec;
		}

		@Override
		public VCFHeader header() { return header; }

		@Override
		public Callable<List<VariantContext>> nextBatch(final int batchSize) throws IOException {
			final ByteArrayOutputStream batch = new ByteArrayOutputStream();
			int count = 0;
			while (count < batchSize) {
				final int n = in.readNBytes(lengths, 0, 8);
				if (n == 0) break; // clean end of file
				if (n < 8) throw new EOFException("Truncated BCF record length fields");
				final ByteBuffer bb = ByteBuffer.wrap(lengths).order(ByteOrder.LITTLE_ENDIAN);
				final long recordLength = Integer.toUnsignedLong(bb.getInt(0)) + Integer.toUnsignedLong(bb.getInt(4));
				if (recordLength > Integer.MAX_VALUE) throw new IOException("BCF record too large: " + recordLength + " bytes");
				final byte[] body = in.readNBytes((int) recordLength);
				if (body.length < recordLength) throw new EOFException("Truncated BCF record");
				batch.write(lengths);
				batch.write(body);
				count++;
			}
			if (count == 0) return null;
			final byte[] bytes = batch.toByteArray();
			final int numRecords = count;
			return () -> {
				final BCF2Codec codec = codecs.get();
				final PositionalBufferedStream stream = new PositionalBufferedStream(new ByteArrayInputStream(bytes), Math.min(bytes.length, 1 << 16));
				final List<VariantContext> result = new ArrayList<>(numRecords);
				for (int i = 0; i < numRecords; i++)
					result.add(fullyDecodeGenotypes(codec.decode(stream)));
				return result;
			};
		}

		@Override
		public void close() throws IOException { in.close(); }
	}

	private class OrderedIterator implements CloseableIterator<VariantContext> {
		private Iterator<VariantContext> current = Collections.emptyIterator();
		private boolean done = false;

		@Override
		public boolean hasNext() {
			while (!current.hasNext()) {
				if (done || closed) return false;
				final Future<List<VariantContext>> f;
				final List<VariantContext> batch;
				final long t0 = System.nanoTime();
				try {
					f = queue.take();
					if (f == END) {
						done = true;
						consumerWaitNanos.addAndGet(System.nanoTime() - t0);
						return false;
					}
					batch = f.get();
				} catch (InterruptedException e) {
					Thread.currentThread().interrupt();
					throw new RuntimeException("Interrupted while reading VCF", e);
				} catch (ExecutionException e) {
					done = true;
					final Throwable cause = e.getCause();
					if (cause instanceof RuntimeException re) throw re;
					if (cause instanceof Error err) throw err;
					throw new RuntimeException(cause);
				}
				consumerWaitNanos.addAndGet(System.nanoTime() - t0);
				current = batch.iterator();
			}
			return true;
		}

		@Override
		public VariantContext next() {
			if (!hasNext()) throw new NoSuchElementException();
			return current.next();
		}

		@Override
		public void close() {
			ParallelVCFReader.this.close();
		}
	}
}
