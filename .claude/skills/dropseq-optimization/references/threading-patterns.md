# Threading patterns and thread balance

## A. Consumer-side counting (reader scales, calling thread at 100%)
Model: `CountBarcodeSequences`. The calling thread only pulls records from the `SamFileMergeUtil` iterator, groups them (Guava `Iterators.partition`, ~1000) onto a **bounded queue** (2 x consumers); consumers run the same per-record method as the single-thread path. One end marker per consumer. On failure: record it, other consumers skip batches so the queue drains, rethrow on the calling thread.
- Shared state: allow list preloaded into `ConcurrentHashMap<String, LongAdder>` (no measurable cost vs a plain map). Per-worker maps are too memory heavy with millions of keys.
- Split rule (`CountBarcodeSequences.getReaderAndCountingThreads()`): <= 3 threads all read and the calling thread counts; >= 4 one counting thread, rest read. Time follows **reader** threads.

## B. Shard-by-cell workers (per-cell state, order-independent)
Model: `SingleCellRnaSeqMetricsCollector` (`SingleCellRnaSeqMetricsAccumulator`). When per-record work is heavy (e.g. overlap detection) and state is per cell, no sort is needed:
- The calling thread filters (cell set, MQ), then routes each record to worker `hash(cell) % W`, in per-worker batches of 1000 on a bounded queue (4 batches per worker). A cell is only ever touched by one worker, so the per-cell collectors need no locks and no merge. Workers keep draining after a failure so the producer never blocks; the failure is rethrown on the calling thread.
- The accumulator takes a `cell -> per-cell collector` map and calls the child directly (0 workers = inline, so the read loop is the same at any thread count). For Picard collectors, the subclass records children in `makeChildCollector` (called from the superclass constructor, so the map field must have no initializer) and skips the `RG` attribute and header read groups; the distributor only does a read-only `Map.get` then `child.acceptRecord`.
- Check the shared collaborators are safe for concurrent reads (gene `OverlapDetector`, shared sets).
- Before threading a tool that sorts by cell/tag, ask whether the sort is needed. Here dropping it was the largest gain (slice: 9.7 min -> 6.8 min at 1 thread, 3.6 min at 4 threads, output identical).
- Thread split (slice sweep, workers were the bottleneck): 1 reader (stock reader, decode on the calling thread) + N-1 workers beat any 2-reader split at every N. Wall 328 s (1 worker), 198 (2), 155 (3), 118 (7); gain flattens near 4 workers because the calling thread (decode + filter + route) becomes the ceiling. A second decode thread cost CPU without a gain.

## C. Reader/writer balance (rules from experiments)
- The calling thread (merge, split, hand-over, per-record logic) is **additional** to the reader and writer pools; count it when sizing CPUs.
- Both sides matter: 1 reader or 1 writer is clearly worse; the optimum is flat in the middle. Default to roughly balanced: readers = (N+1)/2, writers = N - readers (`SplitBamByCell`). One output file: readers `max(1, N-1)` + one async writer (`CorrectScrnaReadPairs`).
- Production shape matters: benchmark with the real number of inputs/outputs (e.g. 8 and 8), not the 1-2 file convenience data. Output count caps useful writer threads (idea: writers = min(N - readers, outputs)).
- Useful ceiling is typically 2-4 threads (writer-bound) unless there are many outputs/inputs.
- Production is often 1 UGER slot, so thread options help less there. The JVM may size pools from the whole node (consider `-XX:ActiveProcessorCount=<slots>` in `JAVA_OPTS`). The launcher appends `-Xmx4g` after `JAVA_OPTS=-Xmx16g` (last wins, not verified).

## D. Memory
Thread count raises RSS (in-flight batches, decode pool, per-file read-ahead, per-output batches). Measure peak RSS with `/usr/bin/time -l` against the Nextflow `-m` setting (production: 80% of VM memory) and give the user the number; try a smaller `-Xmx` to find the real minimum when the heap was at its cap.
