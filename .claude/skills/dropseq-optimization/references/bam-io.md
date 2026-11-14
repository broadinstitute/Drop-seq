# BAM and text I/O paths

## Parallel BAM reader: `SamFileMergeUtil`
Package `org.broadinstitute.dropseqrna.utils.readiterators` (+ `.parallel` internals: `BamDecodePool`, `ParallelBamReader`, `InterleavingSamRecordIterator`, `ReadMode`).

```java
SamHeaderAndIterator h = SamFileMergeUtil.mergeInputPaths(paths, maintainSort, SamReaderFactory.makeDefault(), readerThreads);
// ... h.iterator ...; h.iterator.close() stops the pool.  Overload with explicit ReadMode {AUTO, STOCK, POOL, ASYNC}.
```
- `numThreads == 1` = stock htsjdk reader, no pool. `ReadMode.AUTO`: pool, except one input with an empty sequence dictionary (unaligned) uses htsjdk async IO.
- One shared pool of `numThreads` daemon threads however many files are open. Records are fully decoded on the pool and are fresh objects bound to the caller's header in `next()`. Nothing opens until `iterator()`. `numThreads` counts pool threads: 2 = one reading + one decoding.
- Scaling: grows with the number of files read at once. One aligned file gains up to ~2 threads (compressed) / ~4 (uncompressed); one unaligned file gains little. Beyond ~4 threads rarely helps one tool.
- `maintainSort=true` reads all files concurrently and merges on the consumer thread (can cap scaling).
- `maintainSort=false` with `numThreads > 1` and several inputs **interleaves** inputs round-robin by 250-record batch, so record order depends on thread count. Use only in order-independent tools. Order-sensitive tools (byte-identical BAM output) use `maintainSort=true` or 1 thread.
- Batch 250 records, read-ahead = `numThreads` batches. Records in flight ~500-4000 is the sweet spot; 8000+ falls off a cliff.
- The `parallel` package must not parse BAM or use htsjdk internals (an htsjdk upgrade is coming); differential tests against htsjdk's own reader guard this.
- Wired tools (copy these): `GatherReadQualityMetrics`, `CountBarcodeSequences`, `CorrectScrnaReadPairs`, `SplitBamByCell`, `BamTagHistogram`, `SelectCellsByNumTranscripts`, `DigitalExpression`, `SingleCellRnaSeqMetricsCollector`. Not wired, reader-bound candidates: `DownsampleBamByTag`, `FilterBamByTag` (when few reads pass). Compute-bound tools (DetectDoublets, ACTS, GDAC) do not gain from the reader alone.
- Caller hazards: header-only callers never close the iterator; records retained after `next()` (never pool/reuse records); callers mutate the returned header (workers use a private copy, so add read groups etc. to the returned header before iterating); the same input may be opened several times; some factories request EAGERLY_DECODE; other factory options (reference, inflater factory) are not applied to parallel reads, only validation stringency.

## Codecs
- Picard installs Intel GKL (de)flater by default. `--USE_JDK_DEFLATER true --USE_JDK_INFLATER true` (misleading names) keep htsjdk 5's **libdeflate**. They also select the codec for non-BAM `.gz` written with `ParallelBgzfOutputStream`.
- libdeflate beats Intel at compression level 5; Intel beats it at level 1; level 0 is equal. The inflater flag helps compressed-input read-bound tools, nothing on level-0 BAMs.
- GATK-bundled modules (htsjdk < 5) have no libdeflate.
- Gotcha: if `java.io.tmpdir` is `noexec` (UGER nodes), GKL and jlibdeflate silently fall back to `java.util.zip`. Set `-Djava.io.tmpdir=<exec dir>` and check the banner (`Deflater: Intel; Inflater: Intel`).

## Writing BAMs
Writing a modified record costs 3-6x more than reading it (encode + deflate), so the single writer is usually the ceiling. Options, in order of preference:
1. **Batched writer** (built): `utils/batchwriter/BatchSinkGroup` + `BatchSink<T>` via `utils/SamWriterGroup` (`direct()` or `batched()`; `addOutput`, `write`, `getReadCount`, `close`). A shared pool writes batches; batches of one sink are written in order, different sinks concurrently; max 4 outstanding batches per sink; failures rethrow on the next hand-over/close. Hidden `WRITER_BATCH_SIZE` (default 1000). Parallelism is capped by the number of outputs (one output: more than ~2 threads adds nothing). Records are kept by reference: never mutate after `write`. Memory = outputs x 4 x batch records. Wired into `AbstractSplitBamClp` (`SplitBamByCell`, `CorrectAndSplitScrnaReadPairs`).
2. **htsjdk async writer** (`setUseAsyncIo(true)`): one writer thread, preserves order and bytes; the choice for a single ordered output (`CorrectScrnaReadPairs`). Thrashes with many outputs (one record per hand-over).
3. **Parallel BGZF BAM writer**: designed (`~/bam_parallel_io_plan.md` section 3), not built; the user ruled out custom BGZF writers for now. Reconsider only after measuring with real production codecs on x86. Sorting writers keep a `SortingCollection` in front.

## Order-preserving record-level parallelism: `OrderedConcurrentMapper` + `BatchingIterable`
Model: `TagReadWithGeneFunction.tagReadsConcurrently`. The calling thread reads, `BatchingIterable<T>` groups records, `OrderedConcurrentMapper` maps batches on `numThreads-1` workers, and a consumer thread emits in input order (writer on the consumer).
- Map **batches (~1000)**, never single records. Set `withMaxItems` to a few x workers (default 500k is far too large).
- `--NUM_THREADS N` there = N workers (+ reader + writer). N=1 is the original serial code and the reference for equivalence tests.
- Per-record logic is a small class owned per worker (`withPerThreadFactory`): immutable config shared, mutable metrics per worker and merged at the end (long counts and doubles of whole numbers merge exactly).
- Audit shared singletons for lazy unsynchronized init (`AnnotationUtils` now uses the holder idiom).
- The ceiling is usually the single writer: ~2 workers saturate it.

## Big gzip text reports: `ParallelBgzfOutputStream` (`utils.io`)
```java
OutputStream out = ParallelBgzfOutputStream.openFileForWriting(file, numThreads);   // .gz name -> parallel BGZF, otherwise plain file
// or new ParallelBgzfOutputStream(underlying, level, numThreads[, blockSize])
```
- Output is BGZF (valid gzip): one `BlockCompressedOutputStream` per 1 MiB chunk plus one end-of-file block. Decompressed data equals `GZIPOutputStream`'s; compressed bytes differ.
- htsjdk public API only; no gzip framing, CRC or deflate code, no direct jlibdeflate import (runtime-scope in htsjdk; a direct variant was built and rejected).
- Trap: `BlockCompressedOutputStream.close()` after `close(boolean)` writes a SECOND EOF block; close each chunk exactly once (no try-with-resources).
- `numThreads == 1` compresses inline; more use a pool and the caller writes blocks in order (at most 2 x numThreads in flight). Worker failures are rethrown from `write`/`close`; `close()` is required.
- Wiring: add `--NUM_THREADS`, open the output with `openFileForWriting`, keep a 1 MiB `BufferedOutputStream` between any `PrintStream` and it (an unbuffered `PrintStream` makes one syscall or JNI `deflate` per line), log progress with htsjdk `ProgressLogger` (`record(cellBarcode, 0)` is the convention for non-BAM records). The `--NUM_THREADS` help text and class Javadoc must say the codec flags apply to this output too.
- `IOUtil.slurpLines(File)` does NOT decompress; use `IOUtil.openFileForReading` + `BufferedReader` to read a `.gz` in tests.
